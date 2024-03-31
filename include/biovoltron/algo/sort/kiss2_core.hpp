#pragma once

// TODO: spdlog
// TODO: spdlog print message only when verbose mode is on
// TODO: find out the required headers for the following code
// TODO: add number of threads!!!
// TODO: may produce incorrect result in rare use cases
#include <biovoltron/algo/sort/constant.hpp>
#include <biovoltron/algo/sort/kiss_common.hpp>
#include <biovoltron/algo/sort/structs.hpp>
#include <biovoltron/algo/sort/utils.hpp>
#include <spdlog/spdlog.h>
#include <spdlog/stopwatch.h>

#include <omp.h>
#include <bit>
#include <execution>
#include <span>
#include <thread>
#include <numeric>
#include <ranges>

namespace biovoltron {

namespace kiss {

template<typename char_type, typename size_type>
void encode_S(
  std::ranges::random_access_range auto &S,
  std::ranges::random_access_range auto &lms,
  std::ranges::random_access_range auto &rank,
  std::ranges::random_access_range auto &is_lms_start,
  size_type& tail_reverse_size,
  size_type bits_for_character,
  size_type l,
  size_type m,
  const size_t num_threads,
  vector<ThreadState<size_type>> &states
) {
  auto l_bits = std::bit_width(l - 1);

  vector<size_type> segment_start(num_threads + 1);

#pragma omp parallel num_threads(num_threads)
  {
    auto tid = omp_get_thread_num();
    auto &state = states[tid];
    distribute_workload<size_type>(m, state);

    state.m = 0;

    for (auto i = state.beg; i + 1 < state.beg + state.len; i++) {
      state.m += (lms[i + 1] - lms[i] + l - 1) >> l_bits;
    }

    if (state.len > 0 && state.beg + state.len < m) {
      state.m += (lms[state.beg + state.len] - lms[(state.beg + state.len) - 1] + l - 1) >> l_bits;
    }
  }

  for (auto i = size_t{}; i < num_threads; i++) {
    segment_start[i] = states[i].m;
  }
  std::exclusive_scan(std::begin(segment_start), std::end(segment_start),
                      std::begin(segment_start), size_type{});
  
  size_type n = S.size();
  tail_reverse_size = size_type{};

#pragma omp parallel num_threads(num_threads)
  {
    auto tid = omp_get_thread_num();
    auto &state = states[tid];
    distribute_workload<size_type>(m, state);

    state.m = segment_start[tid];

    for (auto i = state.beg; i < state.beg + state.len; i++) {
      if (i + 1 >= m)
        continue;
      
      auto l_border = lms[i];
      auto r_border = l_border + l;
      for (auto j = size_type{};; j++) {
        if (l_border >= lms[i + 1])
          break;
        
        if (l_border > n - l) {
#pragma omp atomic
          tail_reverse_size++;
        }

        auto encoded_char = size_type{};
        for (auto k = l_border; k < r_border; k++) {
          auto character = (k < n ? S[k] : 0);
          encoded_char = (encoded_char << bits_for_character) | character;
        }
        
        rank[state.m] = encoded_char;
        is_lms_start[state.m] = (j == 0);
        l_border += l;
        r_border += l;
        state.m++;
      }
    }
  }
}

template<typename size_type>
auto get_key(
  std::ranges::random_access_range auto& arr,
  size_type index,
  size_type index_offset
) {
  return (uint64_t(index) + index_offset >= arr.size()) ?
          size_type{} :
          arr[index + index_offset];
}

template<typename size_type>
void initialize_sa(
  std::ranges::random_access_range auto &sa,
  std::ranges::random_access_range auto &rank,
  std::ranges::random_access_range auto &buf,
  size_type tail_reverse_size,
  size_type m_encoded,
  const size_t num_threads,
  vector<ThreadState<size_type>> &states
) {
#pragma omp parallel for num_threads(num_threads)
  for (auto i = size_type{}; i < m_encoded - 1; i++) { sa[i] = i; }

// Since we substitute $ (the sentinel) with infinite amount of the least character,
// the ordering of last few suffixes may be incorrect.
// The swaps ensure that after the radix sort (stable sort),
// the last few suffixes will be moved to the beginning of each segment in reverse order.
  for (auto i = size_type{}; i < tail_reverse_size; i++) {
    if ((m_encoded - 1) - i - 1 > i)
      std::swap(sa[i], sa[(m_encoded - 1) - i - 1]);
  }

  auto bits = sizeof(size_type) * 8;
  auto passes = bits / 16;

  for (auto k = size_type{}; k < passes; k++) {
    auto right_shift_offset = k * 16;

#pragma omp parallel num_threads(num_threads)
    {
      auto tid = omp_get_thread_num();
      auto &state = states[tid];
      distribute_workload<size_type>(m_encoded - 1, state);

      for (auto i = size_type{}; i < KISS2_RADIX_SORT_BUCKET_SIZE; i++)
        state.buffer[i] = size_type{};

      for (auto i = state.beg; i < state.beg + state.len; i++) {
        auto& elem = ((k & 1) ? buf[i] : sa[i]);
        auto key = (rank[elem] >> right_shift_offset) & KISS2_RADIX_SORT_MASK;
        state.buffer[key]++;
      }
    }   

  vector<size_type> segment_start(KISS2_RADIX_SORT_BUCKET_SIZE);
#pragma omp parallel for num_threads(num_threads)
    for (auto i = 0; i < KISS2_RADIX_SORT_BUCKET_SIZE; i++) {
      auto& x = segment_start[i];
      for (auto tid = 0; tid < num_threads; tid++) {
        x += states[tid].buffer[i];
        states[tid].buffer[i] = x - states[tid].buffer[i];
      }
    }

    std::exclusive_scan(std::begin(segment_start), std::end(segment_start),
                        std::begin(segment_start), size_type{});
  
#pragma omp parallel num_threads(num_threads)
    {
      auto tid = omp_get_thread_num();
      auto &state = states[tid];
      distribute_workload<size_type>(m_encoded - 1, state);

      for (auto i = state.beg; i < state.beg + state.len; i++) {
        auto& elem = ((k & 1) ? buf[i] : sa[i]);
        auto key = (rank[elem] >> right_shift_offset) & KISS2_RADIX_SORT_MASK;
        auto index = segment_start[key] + (state.buffer[key]++);
        auto& output_addr = ((k & 1) ? sa[index] : buf[index]);
        output_addr = elem;
      }
    }
  }

  if (passes & 1) {
    std::memmove(sa.data(), buf.data(), (uint64_t)m_encoded * sizeof(size_type));
  }
}

// TODO: check if the extra swap is really needed
template<typename size_type>
void initialize_rank(
  std::ranges::random_access_range auto &sa,
  std::ranges::random_access_range auto &rank,
  std::ranges::random_access_range auto &is_head,
  size_type tail_reverse_size,
  size_type m_encoded,
  const size_t num_threads,
  vector<ThreadState<size_type>> &states
) {
// Define "head" to be an index s.t. its rank differs from
// its previous element.
// buffer[0]: where the block has a head
// buffer[1]: the largest index for a head
#pragma omp parallel num_threads(num_threads)
    {
      auto tid = omp_get_thread_num();
      auto &state = states[tid];
      distribute_workload<size_type>(m_encoded - 1, state);

      auto& is_block_contains_head = state.buffer[0] = false;
      auto& block_largest_head = state.buffer[1] = state.beg;

      for (auto i = state.beg; i < state.beg + state.len; i++) {
        is_head[i] = (i == 0 || 
                      rank[sa[i]] != rank[sa[i - 1]] ||
                      sa[i - 1] >= m_encoded - 1 - tail_reverse_size);
        
        is_block_contains_head |= is_head[i];
        block_largest_head = (is_head[i] ? i + 1 : block_largest_head);
      }
    }

    auto last_block_largest_head = states[0].buffer[1];
    for (auto tid = size_t{1}; tid < num_threads; tid++) {
      auto& is_block_contains_head = states[tid].buffer[0];
      auto& block_largest_head = states[tid].buffer[1];

      if (!is_block_contains_head) {
        block_largest_head = last_block_largest_head;
      }

      last_block_largest_head = block_largest_head;
    }

#pragma omp parallel num_threads(num_threads)
    {
      auto tid = omp_get_thread_num();
      auto &state = states[tid];
      distribute_workload<size_type>(m_encoded - 1, state);

      auto block_last_head = (tid ? states[tid - 1].buffer[1] : size_type{});

      for (auto i = state.beg; i < state.beg + state.len; i++) {
        block_last_head = (is_head[i] ? i + 1 : block_last_head);
        rank[sa[i]] = block_last_head;
      }
    }
}

template<typename size_type>
void sort_sa_blocks(
  std::ranges::random_access_range auto& sa,
  std::ranges::random_access_range auto& rank,
  std::ranges::random_access_range auto& is_head,
  size_type index_offset,
  size_type m_remaining,
  const size_t num_threads,
  vector<ThreadState<size_type>> &states
) {
    auto comparator = [&rank, index_offset](size_type i, size_type j) {
      return get_key(rank, i, index_offset) < get_key(rank, j, index_offset);
    };
#pragma omp parallel num_threads(num_threads)
    {
      auto tid = omp_get_thread_num();
      auto &state = states[tid];
      distribute_workload<size_type>(m_remaining, state);

      auto& last_l_border = state.buffer[0] = state.beg;

      for (auto i = state.beg; i < state.beg + state.len; i++) {
        if (is_head[i]) {
          std::sort(std::begin(sa) + last_l_border, 
                    std::begin(sa) + i, 
                    comparator);
          last_l_border = i;
        }
      }

      std::sort(std::begin(sa) + last_l_border, 
                std::begin(sa) + (state.beg + state.len),
                comparator);

    }
}

template<typename size_type>
void merge_sa_blocks_recursive(
  std::ranges::random_access_range auto& sa,
  std::ranges::random_access_range auto& rank,
  size_type index_offset,
  size_type m_remaining,
  const size_t block_l_border,
  const size_t block_r_border,
  vector<ThreadState<size_type>> &states
) {
  if (block_l_border + 1 == block_r_border)
    return;

  auto block_mid = (block_l_border + block_r_border) / 2;
  merge_sa_blocks_recursive(sa, rank, index_offset, m_remaining, 
                            block_l_border, block_mid, states);
  merge_sa_blocks_recursive(sa, rank, index_offset, m_remaining, 
                            block_mid, block_r_border, states);

  auto& state_l = states[block_l_border];
  auto& l_is_block_contains_head = state_l.buffer[0];
  auto& l_first_segment_end = state_l.buffer[1];
  auto& l_last_segment_head = state_l.buffer[2];

  auto& state_r = states[block_mid];
  auto& r_is_block_contains_head = state_r.buffer[0];
  auto& r_first_segment_end = state_r.buffer[1];
  auto& r_last_segment_head = state_r.buffer[2];
  auto& r_block_start = state_r.beg;

  bool is_indicies_in_range = (l_last_segment_head < m_remaining &&
                               r_block_start < m_remaining &&
                               r_first_segment_end <= m_remaining);
  bool is_ranges_in_order = (l_last_segment_head < r_block_start &&
                             r_block_start < r_first_segment_end);

  if (is_indicies_in_range && is_ranges_in_order) {
    auto comparator = [&rank, index_offset](size_type i, size_type j) {
      return get_key(rank, i, index_offset) < get_key(rank, j, index_offset);
    };

    std::inplace_merge(std::execution::par_unseq,
                       std::begin(sa) + l_last_segment_head,
                       std::begin(sa) + r_block_start,
                       std::begin(sa) + r_first_segment_end,
                       comparator);
  }

  if (!l_is_block_contains_head) {
    l_first_segment_end = r_first_segment_end;
  }
  if (r_is_block_contains_head) {
    l_last_segment_head = r_last_segment_head;
  }
  l_is_block_contains_head |= r_is_block_contains_head;
}

template<typename size_type>
void merge_sa_blocks(
  std::ranges::random_access_range auto& sa,
  std::ranges::random_access_range auto& rank,
  std::ranges::random_access_range auto& is_head,
  size_type index_offset,
  size_type m_remaining,
  const size_t num_threads,
  vector<ThreadState<size_type>> &states
) {
#pragma omp parallel num_threads(num_threads)
    {
      auto tid = omp_get_thread_num();
      auto &state = states[tid];
      distribute_workload<size_type>(m_remaining, state);

      auto &is_block_contains_head = state.buffer[0] = false;
      auto &first_segment_end = state.buffer[1] = state.beg;
      auto &last_segment_head = state.buffer[2] = state.beg;

      for (auto i = state.beg; i < state.beg + state.len; i++) {
        first_segment_end = (is_block_contains_head ? first_segment_end : i);
        is_block_contains_head = (is_block_contains_head || is_head[i]);
        last_segment_head = (is_head[i] ? i : last_segment_head);
      }

      if (!is_block_contains_head && state.beg < m_remaining)
        first_segment_end++;
    }

  merge_sa_blocks_recursive(sa, rank, index_offset, m_remaining, 0, num_threads, states);
}

template<typename size_type>
void
calculate_new_rank_head(
  std::ranges::random_access_range auto& sa,
  std::ranges::random_access_range auto& rank,
  std::ranges::random_access_range auto& is_head,
  size_type index_offset,
  size_type m_remaining,
  const size_t num_threads,
  vector<ThreadState<size_type>> &states
) {
  auto is_new_head = TypeVector(m_remaining, 0);

#pragma omp parallel num_threads(num_threads)
    {
      auto tid = omp_get_thread_num();
      auto &state = states[tid];
      distribute_workload<size_type>(m_remaining, state);

      auto &head_tag = state.buffer[0] = size_type{};
      auto &head_sa_index = state.buffer[1] = size_type{};
      auto &head_mask = state.buffer[2] = size_type{};

      for (auto i = state.beg; i < state.beg + state.len; i++) {
        if (is_head[i]) {
          head_mask |= KISS2_HEAD_MASK;
          head_tag = rank[sa[i]];
          head_sa_index = i;
        } else if (get_key(rank, sa[i], index_offset) !=
                   get_key(rank, sa[i - 1], index_offset)) {
          is_new_head[i] = 1;
          head_mask |= KISS2_NEW_HEAD_MASK;
          head_tag += (i - head_sa_index);
          head_sa_index = i;
        }
      }
    }
 
    auto last_head_tag = states[0].buffer[0];
    auto last_head_sa_index = states[0].buffer[1];
    for (auto i = size_t{1}; i < num_threads; i++) {
      auto& head_tag = states[i].buffer[0];
      auto& head_sa_index = states[i].buffer[1];
      auto head_mask = states[i].buffer[2];

      if (!(head_mask & KISS2_HEAD_MASK)) {
        if (!(head_mask & KISS2_NEW_HEAD_MASK)) {
          head_tag = last_head_tag;
          head_sa_index = last_head_sa_index;
        } else {
          head_tag = last_head_tag + (head_sa_index - last_head_sa_index);
        }
      }

      last_head_tag = head_tag;
      last_head_sa_index = head_sa_index;
    }

#pragma omp parallel num_threads(num_threads)
    {
      auto tid = omp_get_thread_num();
      auto &state = states[tid];
      distribute_workload<size_type>(m_remaining, state);

      auto head_tag = (tid ? states[tid - 1].buffer[0] : size_type{});
      auto head_sa_index = (tid ? states[tid - 1].buffer[1] : size_type{});

      for (auto i = state.beg; i < state.beg + state.len; i++) {
        if (is_head[i]) {
          head_tag = rank[sa[i]];
          head_sa_index = i;
          is_head[i] = 1;
        } else if (is_new_head[i]) {
          auto new_tag = head_tag + (i - head_sa_index);
          rank[sa[i]] = new_tag;
          head_tag = new_tag;
          head_sa_index = i;
          is_head[i] = 1;
        } else {
          rank[sa[i]] = head_tag;
          is_head[i] = 0;
        }
      }
    }
}

template<typename size_type>
size_type
compact(
  std::ranges::random_access_range auto& sa,
  std::ranges::random_access_range auto& buf,
  std::ranges::random_access_range auto& is_head,
  size_type m_remaining,
  const size_t num_threads,
  vector<ThreadState<size_type>> &states
) {
    vector<size_type> compact_start(num_threads + 1);

#pragma omp parallel num_threads(num_threads)
    {
      auto tid = omp_get_thread_num();
      auto &state = states[tid];
      distribute_workload<size_type>(m_remaining, state);

      auto &count_compact = state.buffer[0] = size_type{};

      for (auto i = state.beg; i < state.beg + state.len; i++) {
        bool singleton = is_head[i] && 
                         (i + 1 >= m_remaining || is_head[i + 1]);

        count_compact += (!singleton);
      }
    }

    for (auto i = size_t{}; i < num_threads; i++) {
      compact_start[i] = states[i].buffer[0];
    }

    std::exclusive_scan(std::begin(compact_start), std::end(compact_start),
                        std::begin(compact_start), size_type{});

    auto total_count_compact = compact_start[num_threads];
    auto is_head_buf = TypeVector(total_count_compact, 0);

#pragma omp parallel num_threads(num_threads)
    {
      auto tid = omp_get_thread_num();
      auto &state = states[tid];
      distribute_workload<size_type>(m_remaining, state);

      auto count_compact = compact_start[tid];
      auto count_compact_l_border = compact_start[tid];
      auto count_compact_r_border = compact_start[tid + 1];

      for (auto i = state.beg; i < state.beg + state.len; i++) {
        bool singleton = is_head[i] && 
                         (i + 1 >= m_remaining || is_head[i + 1]);
        if (singleton)
          continue;

        buf[count_compact] = sa[i];
        // prevent race condition in is_head_buf
        if (count_compact - count_compact_l_border < TYPEVECTOR_ELEMENT_SIZE || 
            count_compact_r_border - count_compact < TYPEVECTOR_ELEMENT_SIZE) {
#pragma omp critical
          { is_head_buf[count_compact] = is_head[i]; }
        } else {
          is_head_buf[count_compact] = is_head[i];
        }

        count_compact++;
      }
    }

    std::memmove(sa.data(), buf.data(), (uint64_t)total_count_compact * sizeof(size_type));
    is_head = is_head_buf;

    return total_count_compact;
}

template<typename size_type>
void get_overall_rank(
  std::ranges::random_access_range auto& sa,
  std::ranges::random_access_range auto& rank,
  std::ranges::random_access_range auto& is_head,
  size_type m_remaining,
  const size_t num_threads,
  vector<ThreadState<size_type>> &states
) {

#pragma omp parallel num_threads(num_threads)
    {
      auto tid = omp_get_thread_num();
      auto &state = states[tid];
      distribute_workload<size_type>(m_remaining, state);

      auto& is_block_contains_head = state.buffer[0] = false;
      auto& head_tag = state.buffer[1] = size_type{};
      auto& head_sa_index = state.buffer[2] = size_type{};

      for (auto i = state.beg; i < state.beg + state.len; i++) {
        if (!is_head[i])
          continue;

        is_block_contains_head = true;
        head_tag = rank[sa[i]];
        head_sa_index = i;
      }
    }

    auto last_head_tag = states[0].buffer[1];
    auto last_head_sa_index = states[0].buffer[2];
    for (auto tid = size_t{1}; tid < num_threads; tid++) {
      auto is_block_contains_head = states[tid].buffer[0];
      auto& head_tag = states[tid].buffer[1];
      auto& head_sa_index = states[tid].buffer[2];

      if (!is_block_contains_head) {
        head_tag = last_head_tag;
        head_sa_index = last_head_sa_index;
      }

      last_head_tag = head_tag;
      last_head_sa_index = head_sa_index;
    }

#pragma omp parallel num_threads(num_threads)
    {
      auto tid = omp_get_thread_num();
      auto &state = states[tid];
      distribute_workload<size_type>(m_remaining, state);

      auto head_tag = (tid ? states[tid - 1].buffer[1] : size_type{});
      auto head_sa_index = (tid ? states[tid - 1].buffer[2] : size_type{});

      for (auto i = state.beg; i < state.beg + state.len; i++) {
        if (is_head[i]) {
          head_tag = rank[sa[i]];
          head_sa_index = i;
        } else {
          rank[sa[i]] = head_tag + (i - head_sa_index);
        }
      }
    }
}

template<typename size_type>
void get_overall_sa(
  std::ranges::random_access_range auto& sa,
  std::ranges::random_access_range auto& rank,
  size_type m_encoded,
  const size_t num_threads,
  vector<ThreadState<size_type>> &states
) {
#pragma omp parallel num_threads(num_threads)
    {
      auto tid = omp_get_thread_num();
      auto &state = states[tid];
      distribute_workload<size_type>(m_encoded - 1, state);

      for (auto i = state.beg; i < state.beg + state.len; i++) {
        sa[rank[i] - 1] = i;
      }
    }
}

template<typename size_type>
void place_back_lms(
    std::ranges::random_access_range auto& S,
    std::ranges::random_access_range auto& sa,
    std::ranges::random_access_range auto& buf,
    std::ranges::random_access_range auto& buf2,
    std::ranges::random_access_range auto& is_lms_start,
    const size_type l,
    const size_type m_encoded,
    const size_t num_threads,
    vector<ThreadState<size_type>> &states
  ) {
    auto m = get_lms(S, buf, states);

    auto l_bits = std::bit_width(l - 1);

    vector<size_type> segment_start(num_threads + 1);

#pragma omp parallel num_threads(num_threads)
  {
    auto tid = omp_get_thread_num();
    auto &state = states[tid];
    distribute_workload<size_type>(m, state);

    state.m = 0;

    for (auto i = state.beg; i + 1 < state.beg + state.len; i++) {
      state.m += (buf[i + 1] - buf[i] + l - 1) >> l_bits;
    }

    if (state.len > 0 && state.beg + state.len < m) {
      state.m += (buf[state.beg + state.len] - buf[(state.beg + state.len) - 1] + l - 1) >> l_bits;
    }
  }

  for (auto i = size_t{}; i < num_threads; i++) {
    segment_start[i] = states[i].m;
  }
  std::exclusive_scan(std::begin(segment_start), std::end(segment_start),
                      std::begin(segment_start), size_type{});

#pragma omp parallel num_threads(num_threads)
  {
    auto tid = omp_get_thread_num();
    auto &state = states[tid];
    distribute_workload<size_type>(m, state);

    state.m = segment_start[tid];

    for (auto i = state.beg; i < state.beg + state.len; i++) {
      if (i + 1 >= m)
        continue;
      
      auto l_border = buf[i];
      for (auto j = size_type{};; j++) {
        if (l_border >= buf[i + 1])
          break;

        buf2[state.m++] = l_border;
        l_border += l;
      }
    }
  }

  vector<size_type> count_sorted_lms(num_threads + 1);

#pragma omp parallel num_threads(num_threads)
  {
    auto tid = omp_get_thread_num();
    auto &state = states[tid];
    distribute_workload<size_type>(m_encoded - 1, state);

    state.m = 0;

    for (auto i = state.beg; i < state.beg + state.len; i++) {
      state.m += (is_lms_start[sa[i]]);
    }
  }

  for (auto i = size_t{}; i < num_threads; i++) {
    count_sorted_lms[i] = states[i].m;
  }
  std::exclusive_scan(std::begin(count_sorted_lms), std::end(count_sorted_lms),
                      std::begin(count_sorted_lms), size_type{});

  #pragma omp parallel num_threads(num_threads)
  {
    auto tid = omp_get_thread_num();
    auto &state = states[tid];
    distribute_workload<size_type>(m_encoded - 1, state);

    state.m = count_sorted_lms[tid];

    for (auto i = state.beg; i < state.beg + state.len; i++) {
      if (is_lms_start[sa[i]])
          buf[state.m++] = buf2[sa[i]];
    }
  }

  std::memmove(sa.data() + 1, buf.data(), (uint64_t)(m - 1) * sizeof(size_type));
  auto n = S.size();
  sa[0] = n;
}

template<typename char_type, typename size_type>
void prefix_doubling(
  std::ranges::random_access_range auto &S,
  std::ranges::random_access_range auto &SA,
  size_type k,
  size_type bits_for_character,
  size_type l,
  size_type m,
  size_type m_encoded,
  const size_t num_threads,
  vector<ThreadState<size_type>> &states
) {
  auto sw = spdlog::stopwatch{};
  auto lms = std::ranges::subrange(std::begin(SA),
                                   std::begin(SA) + uint64_t(m_encoded));
  auto rank = std::ranges::subrange(std::begin(SA) + uint64_t(m_encoded),
                                    std::begin(SA) + uint64_t(m_encoded) * 2);
  auto buffer = std::ranges::subrange(std::begin(SA) + uint64_t(m_encoded) * 2,
                                      std::begin(SA) + uint64_t(m_encoded) * 3);
  auto is_lms_start = TypeVector(m_encoded, 0);
  auto tail_reverse_size = size_type{};
  kiss::encode_S<size_type>(S, lms, rank, is_lms_start, tail_reverse_size, bits_for_character, l, m, num_threads, states);
  SPDLOG_DEBUG("encode_S elapsed {}", sw);

  sw = spdlog::stopwatch{};
  auto sa = lms;
  initialize_sa<size_type>(sa, rank, buffer, tail_reverse_size, m_encoded, num_threads, states);
  SPDLOG_DEBUG("initialize_sa elapsed {}", sw);

  sw = spdlog::stopwatch{};
// This Boolean array is used for easier implementation
  auto is_head = TypeVector(m_encoded, 0);
  initialize_rank<size_type>(sa, rank, is_head, tail_reverse_size, m_encoded, num_threads, states);
  SPDLOG_DEBUG("initialize_rank elapsed {}", sw);

  auto m_remaining = m_encoded - 1;
  // TODO: index_offset may overflow
  for (size_type index_offset = 1; index_offset < k; index_offset *= 2) {
    sw = spdlog::stopwatch{};
    sort_sa_blocks(sa, rank, is_head, index_offset, m_remaining, num_threads, states);
    SPDLOG_DEBUG("sort_sa_blocks elapsed {}", sw);

    sw = spdlog::stopwatch{};
    merge_sa_blocks(sa, rank, is_head, index_offset, m_remaining, num_threads, states);
    SPDLOG_DEBUG("merge_sa_blocks elapsed {}", sw);

    sw = spdlog::stopwatch{};
    calculate_new_rank_head(sa, rank, is_head, index_offset, m_remaining, num_threads, states);
    SPDLOG_DEBUG("calculate_new_rank_head elapsed {}", sw);

    sw = spdlog::stopwatch{};
    m_remaining = compact<size_type>(sa, buffer, is_head, m_remaining, num_threads, states);
    SPDLOG_DEBUG("compact elapsed {}", sw);

    if (m_remaining <= 1)
      break;
  }

  sw = spdlog::stopwatch{};
  get_overall_rank<size_type>(sa, rank, is_head, m_remaining, num_threads, states);
  SPDLOG_DEBUG("get_overall_rank elapsed {}", sw);

  sw = spdlog::stopwatch{};
  get_overall_sa<size_type>(sa, rank, m_encoded, num_threads, states);
  SPDLOG_DEBUG("get_overall_sa elapsed {}", sw);

  sw = spdlog::stopwatch{};
  place_back_lms<size_type>(S, sa, rank, buffer, is_lms_start, l, m_encoded, num_threads, states);
  SPDLOG_DEBUG("place_back_lms elapsed {}", sw);
}

template<typename size_type>
auto get_encoded_S_length(
  std::ranges::random_access_range auto& lms,
  const size_type m,
  const size_type l,
  const size_t num_threads,
  vector<ThreadState<size_type>> &states
) {
  auto l_bits = std::bit_width(l - 1);

#pragma omp parallel num_threads(num_threads)
  {
    auto tid = omp_get_thread_num();
    auto &state = states[tid];
    distribute_workload<size_type>(m, state);

    state.m = 0;

    for (auto i = state.beg; i + 1 < state.beg + state.len; i++) {
      state.m += (lms[i + 1] - lms[i] + l - 1) >> l_bits;
    }

    if (state.len > 0 && state.beg + state.len < m) {
      state.m += (lms[state.beg + state.len] - lms[(state.beg + state.len) - 1] + l - 1) >> l_bits;
    }
  }

  auto length = size_type{};
  for (auto &state : states)
    length += state.m;

  length++; // add the length of the sentinel
  return length;
}

template <typename char_type, typename size_type>
void kiss2_suffix_array_internal(
  vector<char_type> &S,
  vector<size_type> &SA,
  size_type k = 256u,
  const size_type character_size = CHAR_SIZE,
  const size_t num_threads = std::thread::hardware_concurrency()
) {
  auto sw = spdlog::stopwatch{};
  if (S.size() == 0)
    return SA = {0}, void();

  size_type n = S.size();
  SA.reserve(n + 1);

  auto states = vector<ThreadState<size_type>>(num_threads);
  SPDLOG_DEBUG("Preparation elapsed {}", sw);

  sw = spdlog::stopwatch{};
  SA.resize(n + 1);
  SPDLOG_DEBUG("SA.resize(n + 1) elapsed {}", sw);

  sw = spdlog::stopwatch{};
  auto m = get_lms(S, SA, states);
  SPDLOG_DEBUG("get_lms elapsed {}", sw);

  sw = spdlog::stopwatch{};
  const size_type bits_for_character = std::bit_width(character_size - 1);
  const size_type l = 8 * sizeof(size_type) / bits_for_character;
  auto m_encoded = get_encoded_S_length(SA, m, l, num_threads, states);
  SPDLOG_DEBUG("get_encoded_S_length elapsed {}", sw);

  if (uint64_t(m_encoded) * 3 > n + 1) {
    sw = spdlog::stopwatch{};
    SA.resize(uint64_t(m_encoded) * 3);
    SPDLOG_DEBUG("SA.resize(3 * m_encoded) elapsed {}", sw);
  }

  sw = spdlog::stopwatch{};
  prefix_doubling<char_type, size_type>(S, SA, k, bits_for_character, l, m, m_encoded, num_threads, states);
  SPDLOG_DEBUG("prefix_doubling elapsed {}", sw);

  sw = spdlog::stopwatch{};
  SA.resize(n + 1);
  auto SA1 = std::ranges::subrange(SA.begin() + 1, SA.begin() + m);
  put_lms_suffix(S, SA, SA1, states);
  SPDLOG_DEBUG("put_lms_suffix elapsed {}", sw);

  sw = spdlog::stopwatch{};
  induced_sort(S, SA, states);
  SPDLOG_DEBUG("induced_sort elapsed {}", sw);
}

template <typename char_type, typename size_type>
void kiss2_suffix_array_dna(
  vector<char_type> &S,
  vector<size_type> &SA,
  size_type k = 256u,
  const size_t num_threads = std::thread::hardware_concurrency()
) {
    kiss2_suffix_array_internal<char_type, size_type>(S, SA, k, 4, num_threads);
}

template <typename char_type, typename size_type>
void kiss2_suffix_array(
  vector<char_type> &S,
  vector<size_type> &SA,
  size_type k = 256u,
  const size_t num_threads = std::thread::hardware_concurrency()
) {
    kiss2_suffix_array_internal<char_type, size_type>(S, SA, k, CHAR_SIZE, num_threads);
}

} // namespace kiss

} // namespace biovoltron