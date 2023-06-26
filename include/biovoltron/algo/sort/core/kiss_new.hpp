// FIXME: refactor the code!
#pragma once

#include <bits/stdc++.h>
#include <spdlog/spdlog.h>
#include <spdlog/stopwatch.h>

#include <algorithm>
#include <biovoltron/container/xbit_vector.hpp>
#include <biovoltron/utility/istring.hpp>
#include <biovoltron/utility/thread_pool.hpp>
#include <chrono>
#include <cstdlib>
#include <execution>
#include <numeric>
#include <tuple>
#include <vector>

#include "omp.h"

namespace biovoltron {

namespace kiss {

#define MASK        ((1 << 16) - 1)
#define NUM_THREADS std::thread::hardware_concurrency()
#define PREFIX_DOUBLING_NUM_THREADS(n) \
  std::min((uint32_t)NUM_THREADS, (uint32_t)n)

enum SUFFIX_TYPE { L_TYPE = 0, S_TYPE = 1 };
enum BLOCK_ELEM_TYPE { NONHEAD = 0, HEAD = 1 };
using TypeVector = biovoltron::detail::XbitVector<1, std::uint8_t,
                                                  std::allocator<std::uint8_t>>;

template<typename size_type>
struct MergeBlockData {
  std::vector<size_type> block_have_head, prev_block_end, terminal_block_head;
  MergeBlockData(size_type threads) {
    block_have_head.resize(threads);
    prev_block_end.resize(threads);
    terminal_block_head.resize(threads);
  }
};

template<typename size_type>
auto
is_LMS(const auto& T, size_type i) {
  auto n = (size_type)T.size();
  return i == n
         or (i > 0 and T[i - 1] == SUFFIX_TYPE::L_TYPE
             and T[i] == SUFFIX_TYPE::S_TYPE);
}

auto
get_type_block_range(auto num_items, auto num_blocks, auto block_idx) {
  // in order to make TypeVector(a.k.a XBitVector) thread-safe,
  // the element should be group by 8 bit at least
  auto N = (num_items + 7) / 8, m = num_blocks;
  auto counts = N / m, remain = N % m;
  auto L = counts * block_idx + std::min(block_idx, remain);
  auto R = L + counts + (block_idx < remain);
  return std::tuple{std::min(num_items, L * 8), std::min(num_items, R * 8)};
}

template<typename size_type>
auto
len_compressed_string(const auto& T, size_type compress_block_length) {
  auto n = (size_type)T.size();
  std::vector<size_type> appear(NUM_THREADS), first_place(NUM_THREADS),
    last_place(NUM_THREADS);
  auto len = size_type{};

#pragma omp parallel for num_threads(NUM_THREADS)
  for (auto tid = size_type{}; tid < NUM_THREADS; tid++) {
    auto [L, R] = get_type_block_range(n, NUM_THREADS, tid);
    bool cur_appear = false;
    auto cur_first = R, cur_last = L;
    auto local_len = size_type{};
    for (auto i = L; i < R; i++) {
      if (!is_LMS(T, i))
        continue;
      if (cur_appear) {
        local_len
          += (i - cur_last + compress_block_length) / compress_block_length;
        cur_last = i;
      } else {
        cur_first = i;
        cur_last = i;
      }
      cur_appear = true;
    }
#pragma omp critical
    {
      len += local_len;
      appear[tid] = cur_appear;
      first_place[tid] = cur_first;
      last_place[tid] = cur_last;
    }
  }
  auto last_lms = n;
  for (auto i = size_type{NUM_THREADS - 1}; ~i; i--) {
    if (appear[i]) {
      len += (last_lms - last_place[i] + compress_block_length)
             / compress_block_length;
      last_lms = first_place[i];
    }
  }
  return len;
}

template<typename size_type>
void
get_lms_indices_in_new_string(auto& lms,
                              const std::ranges::random_access_range auto& buf,
                              size_type compress_block_length) {
  auto n1 = (size_type)lms.size() - 1;
  std::vector<size_type> block_start(NUM_THREADS);
  // calculate the entry occupied for each LMS substring
#pragma omp parallel for num_threads(NUM_THREADS)
  for (auto tid = size_type{}; tid < NUM_THREADS; tid++) {
    auto [L, R] = get_type_block_range(n1 + 1, NUM_THREADS, tid);
    for (auto i = L; i < R; i++) {
      buf[i] = (i == n1 ? 1 :
                          (lms[i + 1] - lms[i] + compress_block_length)
                            / compress_block_length);
      if (i > L)
        buf[i] += buf[i - 1];
    }
#pragma omp critical
    { block_start[tid] = buf[R - 1]; }
  }
  // prefix sum
  std::inclusive_scan(std::begin(block_start), std::end(block_start),
                      std::begin(block_start));
// add to buf
#pragma omp parallel for num_threads(NUM_THREADS)
  for (auto tid = size_type{1}; tid < NUM_THREADS; tid++) {
    auto [L, R] = get_type_block_range(n1 + 1, NUM_THREADS, tid);
    for (auto i = L; i < R; i++) { buf[i] += block_start[tid - 1]; }
  }
}

template<typename size_type>
void
compress_string(const std::ranges::random_access_range auto& S, auto& lms,
                auto& rank, const std::ranges::random_access_range auto& buf,
                const std::ranges::random_access_range auto& starting_position,
                std::ranges::random_access_range auto& valid_position,
                size_type K, size_type compress_block_length) {
  auto n = (size_type)S.size();
  auto n2 = (size_type)rank.size() - 1;
  get_lms_indices_in_new_string(lms, buf, compress_block_length);

// calculate the compressed strings
#pragma omp parallel for num_threads(NUM_THREADS)
  for (auto tid = size_type{}; tid < NUM_THREADS; tid++) {
    auto [L, R] = get_type_block_range(n2 + 1, NUM_THREADS, tid);
    // binary search the first index
    size_type idx
      = std::upper_bound(std::begin(buf), std::end(buf), L) - std::begin(buf);
    for (auto i = L; i < R; i++) {
      if (buf[idx] <= i)
        idx++;
      auto offset_count = i - (idx > 0 ? buf[idx - 1] : 0);
      auto l_border = lms[idx] + offset_count * compress_block_length;
      auto r_border = l_border + compress_block_length;

      auto compressed_value = size_type{};
      for (auto j = l_border; j < r_border; j++) {
        // FIXME: dirty fix!
        auto character_value = (j < n ? S[j] + 1 : 0);
        compressed_value = compressed_value * K + character_value;
      }
      rank[i] = compressed_value;
      starting_position[i] = l_border;
      valid_position[i] = (offset_count == 0);
    }
  }
}

template<typename size_type>
void
init_counting_sort(auto& elems_i, auto& elems_o, auto&& member,
                   auto right_shift_offset) {
  auto offsets = std::vector<std::vector<size_type>>(
    NUM_THREADS, std::vector<size_type>(1 << 16, 0));

#pragma omp parallel for num_threads(NUM_THREADS)
  for (auto i = size_type{}; i < elems_i.size(); i++) {
    auto tid = omp_get_thread_num();
    auto& elem = elems_i[i];
    offsets[tid][(member[elem] >> right_shift_offset) & MASK]++;
  }

  // FIXME: dirty fix to make NUM_THREADS constant
  auto total = std::vector<size_type>(1 << 16, 0);
  const size_type num_threads = NUM_THREADS;
  for (auto i = 0; i < (1 << 16); i++) {
    auto& x = total[i];
    for (auto tid = 0; tid < num_threads; tid++) {
      x += offsets[tid][i];
      offsets[tid][i] = x - offsets[tid][i];
    }
  }

  std::exclusive_scan(total.begin(), total.end(), total.begin(), size_type{},
                      std::plus<>{});

#pragma omp parallel for num_threads(NUM_THREADS)
  for (auto tid = 0; tid < NUM_THREADS; tid++) {
    for (auto i = 0; i < (1 << 16); i++) offsets[tid][i] += total[i];
  }

#pragma omp parallel for num_threads(NUM_THREADS)
  for (auto i = size_type{}; i < elems_i.size(); i++) {
    auto tid = omp_get_thread_num();
    auto& elem = elems_i[i];
    auto idx = offsets[tid][(member[elem] >> right_shift_offset) & MASK]++;
    elems_o[idx] = elem;
  }
}

template<typename size_type>
void
init_radix_sort(auto& elems_i, auto& elems_o, auto&& member) {
  init_counting_sort<size_type>(elems_i, elems_o, member, 0);
  init_counting_sort<size_type>(elems_o, elems_i, member, 16);
}

template<typename size_type>
void
init_sa(const std::ranges::random_access_range auto& sa,
        const std::ranges::random_access_range auto& rank,
        const std::ranges::random_access_range auto& buf) {
  auto n2 = (size_type)sa.size() - 1;
#pragma omp parallel for num_threads(NUM_THREADS)
  for (auto i = size_type{}; i < n2 + 1; i++) { sa[i] = i; }
  init_radix_sort<size_type>(sa, buf, rank);
}

template<typename size_type>
void
init_rank(const std::ranges::random_access_range auto& sa,
          std::ranges::random_access_range auto& rank,
          std::ranges::random_access_range auto& buf,
          std::ranges::random_access_range auto& is_head) {
  auto n2 = (size_type)sa.size() - 1;
  std::vector<bool> have_label(NUM_THREADS);
  std::vector<size_type> last_label(NUM_THREADS);
  // pre scan
#pragma omp parallel for num_threads(NUM_THREADS)
  for (auto tid = size_type{}; tid < NUM_THREADS; tid++) {
    auto [L, R] = get_type_block_range(n2 + 1, NUM_THREADS, tid);
    auto cur_last_label = L;
    bool cur_have_label = false;
    for (auto i = L; i < R; i++) {
      if (i == 0 || rank[sa[i]] != rank[sa[i - 1]]) {
        is_head[i] = BLOCK_ELEM_TYPE::HEAD;
        cur_have_label = true;
        cur_last_label = i;
      }
    }
#pragma omp critical
    {
      have_label[tid] = cur_have_label;
      last_label[tid] = cur_last_label;
    }
  }
  // collect last labels
  for (auto tid = size_type{1}; tid < NUM_THREADS; tid++) {
    if (!have_label[tid]) {
      last_label[tid] = last_label[tid - 1];
    }
  }
  // post scan
#pragma omp parallel for num_threads(NUM_THREADS)
  for (auto tid = size_type{}; tid < NUM_THREADS; tid++) {
    auto [L, R] = get_type_block_range(n2 + 1, NUM_THREADS, tid);
    auto cur_last_label = (tid > 0 ? last_label[tid - 1] : 0);
    for (auto i = L; i < R; i++) {
      if (i == 0 || rank[sa[i]] != rank[sa[i - 1]]) {
        buf[sa[i]] = i;
        cur_last_label = i;
      } else {
        buf[sa[i]] = cur_last_label;
      }
    }
  }
  // extra swap
  // #pragma omp parallel for num_threads(NUM_THREADS)
  //   for (auto i = size_type{}; i < n2 + 1; i++) std::swap(buf[i], rank[i]);
  auto new_rank
    = std::ranges::subrange(std::begin(buf), std::begin(buf) + (n2 + 1));
  buf = rank;
  rank = new_rank;
}

template<typename size_type>
auto
get_key(const std::ranges::random_access_range auto& rank, size_type idx,
        size_type index_offset) {
  return ((uint64_t)idx + index_offset >= rank.size()) ?
           size_type{} :
           rank[idx + index_offset];
}

template<typename size_type>
void
radix_sort_bit(const std::ranges::random_access_range auto& sa,
               const std::ranges::random_access_range auto& rank,
               const std::ranges::random_access_range auto& buf,
               size_type index_offset, size_type L, size_type R,
               int bit_offset) {
  std::vector<size_type> start((1 << 16) + 1);
  for (auto i = L; i < R; i++) {
    size_type key = get_key(rank, sa[i], index_offset);
    size_type cur_index = (key >> bit_offset) & MASK;
    start[cur_index + 1]++;
  }
  std::inclusive_scan(std::begin(start), std::end(start), std::begin(start));
  for (auto i = L; i < R; i++) {
    size_type key = get_key(rank, sa[i], index_offset);
    size_type cur_index = (key >> bit_offset) & MASK;
    buf[L + start[cur_index]] = sa[i];
    start[cur_index]++;
  }
  std::swap_ranges(std::begin(sa) + L, std::begin(sa) + R, std::begin(buf) + L);
}

template<typename size_type>
void
sort_same_sa_value(const std::ranges::random_access_range auto& sa,
                   const std::ranges::random_access_range auto& rank,
                   const std::ranges::random_access_range auto& buf,
                   size_type index_offset, size_type L, size_type R) {
  if (R - L <= (1 << 10)) {  // threshold can be adjusted
    std::stable_sort(std::begin(sa) + L, std::begin(sa) + R,
                     [&](size_type idx1, size_type idx2) {
                       return get_key(rank, idx1, index_offset)
                              < get_key(rank, idx2, index_offset);
                     });
  } else {
    radix_sort_bit(sa, rank, buf, index_offset, L, R, 0);
    radix_sort_bit(sa, rank, buf, index_offset, L, R, 16);
  }
}

template<typename size_type>
void
sort_sa_blocks(const std::ranges::random_access_range auto& sa,
               const std::ranges::random_access_range auto& rank,
               const std::ranges::random_access_range auto& buf,
               const std::ranges::random_access_range auto& is_head,
               size_type index_offset) {
#pragma omp parallel for num_threads(NUM_THREADS)
  for (auto tid = size_type{}; tid < NUM_THREADS; tid++) {
    auto n_ = (size_type)sa.size() - 1;
    auto [L, R] = get_type_block_range(n_ + 1, NUM_THREADS, tid);
    auto last_left = L;
    for (auto i = L; i < R; i++) {
      if (is_head[i]) {
        sort_same_sa_value(sa, rank, buf, index_offset, last_left, i);
        last_left = i;
      }
    }
    sort_same_sa_value(sa, rank, buf, index_offset, last_left, R);
  }
}

template<typename size_type>
std::tuple<bool, size_type, size_type>
merge_sa_blocks_recursive(const std::ranges::random_access_range auto& sa,
                          const std::ranges::random_access_range auto& rank,
                          const std::ranges::random_access_range auto& buf,
                          MergeBlockData<size_type>& merge_block_data,
                          size_type index_offset, size_type L, size_type R) {
  if (L + 1 == R) {
    return std::make_tuple(merge_block_data.block_have_head[L],
                           merge_block_data.prev_block_end[L],
                           merge_block_data.terminal_block_head[L]);
  }
  auto mid = (L + R) / 2;
  auto [l_have_head, l_prev_block_end, l_terminal_block_head]
    = merge_sa_blocks_recursive(sa, rank, buf, merge_block_data, index_offset,
                                L, mid);
  auto [r_have_head, r_prev_block_end, r_terminal_block_head]
    = merge_sa_blocks_recursive(sa, rank, buf, merge_block_data, index_offset,
                                mid, R);

  auto n_ = (size_type)sa.size() - 1;
  auto [mid_block_start, _] = get_type_block_range(n_ + 1, NUM_THREADS, mid);
  // FIXME: not sure if it is correct
  if (l_terminal_block_head < n_ + 1 && mid_block_start < n_ + 1
      && r_prev_block_end <= n_ + 1 && l_terminal_block_head < mid_block_start
      && mid_block_start < r_prev_block_end) {
    std::merge(std::execution::par, std::begin(sa) + l_terminal_block_head,
               std::begin(sa) + mid_block_start,
               std::begin(sa) + mid_block_start,
               std::begin(sa) + r_prev_block_end,
               std::begin(buf) + l_terminal_block_head,
               [&](size_type idx1, size_type idx2) {
                 return get_key(rank, idx1, index_offset)
                        < get_key(rank, idx2, index_offset);
               });
    std::swap_ranges(std::begin(sa) + l_terminal_block_head,
                     std::begin(sa) + r_prev_block_end,
                     std::begin(buf) + l_terminal_block_head);
  }

  auto new_have_head = (l_have_head || r_have_head);
  auto new_prev_block_end = (l_have_head ? l_prev_block_end : r_prev_block_end);
  auto new_terminal_block_head
    = (r_have_head ? r_terminal_block_head : l_terminal_block_head);
  return std::make_tuple(new_have_head, new_prev_block_end,
                         new_terminal_block_head);
}

template<typename size_type>
void
merge_sa_blocks(const std::ranges::random_access_range auto& sa,
                const std::ranges::random_access_range auto& rank,
                const std::ranges::random_access_range auto& buf,
                const std::ranges::random_access_range auto& is_head,
                size_type index_offset) {
  MergeBlockData merge_block_data(NUM_THREADS);
#pragma omp parallel for num_threads(NUM_THREADS)
  for (auto tid = size_type{}; tid < NUM_THREADS; tid++) {
    auto n_ = (size_type)sa.size() - 1;
    auto [L, R] = get_type_block_range(n_ + 1, NUM_THREADS, tid);

    size_type local_prev_block_end = L, local_terminal_block_head = L;
    bool have_true = false;
    for (size_type i = L; i < R; i++) {
      if (!have_true)
        local_prev_block_end = i;
      if (is_head[i]) {
        have_true = true;
        local_terminal_block_head = i;
      }
    }
    if (!have_true && L < n_ + 1)
      local_prev_block_end++;
#pragma omp critical
    {
      merge_block_data.block_have_head[tid] = have_true;
      merge_block_data.prev_block_end[tid] = local_prev_block_end;
      merge_block_data.terminal_block_head[tid] = local_terminal_block_head;
    }
  }
  merge_sa_blocks_recursive(sa, rank, buf, merge_block_data, index_offset,
                            size_type{}, NUM_THREADS);
}

template<typename size_type>
void
calculate_new_rank_head(const std::ranges::random_access_range auto& sa,
                        std::ranges::random_access_range auto& rank,
                        std::ranges::random_access_range auto& buf,
                        std::ranges::random_access_range auto& is_head,
                        size_type index_offset) {
  // for (auto i : rank) std::cout << i << ' ';
  // std::cout << std::endl;
  // std::cout << get_key<size_type>(rank, 0, 0) << std::endl;
  // std::cout << get_key<size_type>(rank, 1, 0) << std::endl;
  auto n_ = (size_type)sa.size() - 1;
  auto n2 = (size_type)rank.size() - 1;
  auto is_new_head = TypeVector(n_, SUFFIX_TYPE::L_TYPE);
  // std::copy(std::execution::par, std::begin(rank), std::begin(rank) + (n2 +
  // 1), std::begin(buf));
  std::vector<size_type> head_tag(NUM_THREADS), head_sa_index(NUM_THREADS);
  std::vector<int> determine_level(NUM_THREADS);
#pragma omp parallel for num_threads(NUM_THREADS)
  for (auto tid = size_type{}; tid < NUM_THREADS; tid++) {
    auto [L, R] = get_type_block_range(n_ + 1, NUM_THREADS, tid);
    size_type local_head_tag = size_type{}, local_head_sa_index = size_type{};
    int local_determine_level = 0;
    for (auto i = L; i < R; i++) {
      if (is_head[i]) {
        local_determine_level = 2;
        local_head_tag = rank[sa[i]];
        buf[i] = local_head_tag;
        local_head_sa_index = i;
      } else if (get_key<size_type>(rank, sa[i], index_offset)
                 != get_key<size_type>(rank, sa[i - 1], index_offset)) {
        is_new_head[i] = 1;
        local_determine_level = std::max(local_determine_level, 1);
        local_head_tag = (i - local_head_sa_index) + local_head_tag;
        local_head_sa_index = i;
      }
    }
#pragma omp critical
    {
      head_tag[tid] = local_head_tag;
      head_sa_index[tid] = local_head_sa_index;
      determine_level[tid] = local_determine_level;
    }
  }
  for (auto i = size_type{1}; i < NUM_THREADS; i++) {
    if (determine_level[i] == 0) {
      head_tag[i] = head_tag[i - 1];
      head_sa_index[i] = head_sa_index[i - 1];
    } else if (determine_level[i] == 1) {
      head_tag[i] = head_tag[i - 1] + (head_sa_index[i] - head_sa_index[i - 1]);
    }
  }
#pragma omp parallel for num_threads(NUM_THREADS)
  for (auto tid = size_type{}; tid < NUM_THREADS; tid++) {
    auto [L, R] = get_type_block_range(n_ + 1, NUM_THREADS, tid);
    size_type local_head_tag = (tid ? head_tag[tid - 1] : size_type{});
    size_type local_head_sa_index
      = (tid ? head_sa_index[tid - 1] : size_type{});
    for (auto i = L; i < R; i++) {
      if (is_head[i]) {
        local_head_tag = buf[i];
        local_head_sa_index = i;
        is_head[i] = BLOCK_ELEM_TYPE::HEAD;
      } else if (is_new_head[i]) {
        auto new_tag = local_head_tag + (i - local_head_sa_index);
        rank[sa[i]] = new_tag;
        local_head_tag = new_tag;
        local_head_sa_index = i;
        is_head[i] = BLOCK_ELEM_TYPE::HEAD;
      } else {
        rank[sa[i]] = local_head_tag;
        is_head[i] = BLOCK_ELEM_TYPE::NONHEAD;
      }
    }
  }
  // auto new_rank
  //   = std::ranges::subrange(std::begin(buf), std::begin(buf) + (n2 + 1));
  // buf = rank;
  // rank = new_rank;
}

template<typename size_type>
void
compact(std::ranges::random_access_range auto& sa,
        const std::ranges::random_access_range auto& rank,
        std::ranges::random_access_range auto& buf,
        std::ranges::random_access_range auto& is_head,
        const std::ranges::random_access_range auto& starting_position,
        size_type index_offset, size_type sort_len) {
  auto sw2 = spdlog::stopwatch{};
  auto n_ = (size_type)sa.size() - 1;
  auto n2 = (size_type)rank.size() - 1;
  std::vector<size_type> count_compat(NUM_THREADS + 1);
#pragma omp parallel for num_threads(NUM_THREADS)
  for (auto tid = size_type{}; tid < NUM_THREADS; tid++) {
    size_type count_compat_cur = 0;
    auto [L, R] = get_type_block_range(n_ + 1, NUM_THREADS, tid);

    for (auto i = L; i < R; i++) {
      bool singleton = is_head[i] && (i >= n_ || is_head[i + 1]);
      bool length_exceeded
        = (i + index_offset > n2 + 1)
          || (starting_position[i + index_offset] - starting_position[i]
              >= sort_len);
      if (!singleton && !length_exceeded)
        count_compat_cur++;
    }
#pragma omp critical
    { count_compat[tid + 1] = count_compat_cur; }
  }
  std::inclusive_scan(std::begin(count_compat), std::end(count_compat),
                      std::begin(count_compat));
  SPDLOG_DEBUG("compact 1 elapsed {}", sw2);
  sw2 = spdlog::stopwatch{};

  auto new_n_ = count_compat[NUM_THREADS];
  auto is_head_buf = TypeVector(new_n_, BLOCK_ELEM_TYPE::NONHEAD);
#pragma omp parallel for num_threads(NUM_THREADS)
  for (auto tid = size_type{}; tid < NUM_THREADS; tid++) {
    size_type count_compat_cur = count_compat[tid];
    auto [L, R] = get_type_block_range(n_ + 1, NUM_THREADS, tid);

    for (auto i = L; i < R; i++) {
      bool singleton = is_head[i] && (i >= n_ || is_head[i + 1]);
      bool length_exceeded
        = (i + index_offset > n2 + 1)
          || (starting_position[i + index_offset] - starting_position[i]
              >= sort_len);
      if (!singleton && !length_exceeded) {
        buf[count_compat_cur] = sa[i];               // (!)
        is_head_buf[count_compat_cur] = is_head[i];  // (!)
        count_compat_cur++;
      }
    }
  }
  SPDLOG_DEBUG("compact 2 elapsed {}", sw2);
  // #pragma omp parallel num_threads(NUM_THREADS)
  //   for (auto i = size_type{}; i < new_n_; i++) { sa[i] = buf[i]; }
  auto new_buf
    = std::ranges::subrange(std::begin(sa), std::begin(sa) + (n2 + 1));
  sa = buf;
  buf = new_buf;
  sa = std::ranges::subrange(std::begin(sa), std::begin(sa) + (new_n_));
  is_head = is_head_buf;
}

template<typename size_type>
void
radix_sort(std::ranges::random_access_range auto& sa,
           std::ranges::random_access_range auto& rank,
           std::ranges::random_access_range auto& buf,
           std::ranges::random_access_range auto& is_head,
           const std::ranges::random_access_range auto& starting_position,
           size_type index_offset, size_type sort_len) {
  // std::cout << "<<\n";
  //   for (int i = 0; i < 50; i++) std::cout << rank[i] << ' ';
  //   std::cout << std::endl;
  //   for (int i = 0; i < 50; i++) std::cout << sa[i] << ' ';
  //   std::cout << std::endl;
  //   bool ok = true;
  // #pragma omp parallel for num_threads(NUM_THREADS)
  //   for (auto i = size_type{}; i < sa.size() - 1; i++) {
  //     if (sa[i] == sa[i + 1]) {
  // #pragma omp critical
  //       { ok = false; }
  //     }
  //   }
  //   assert(ok);
  //   std::cout << ok << std::endl;
  auto sw2 = spdlog::stopwatch{};
  sort_sa_blocks(sa, rank, buf, is_head, index_offset);
  SPDLOG_DEBUG("sort sa blocks elapsed {}", sw2);
  // #pragma omp parallel for num_threads(NUM_THREADS)
  //   for (auto i = size_type{}; i < sa.size() - 1; i++) {
  //     if (sa[i] == sa[i + 1]) {
  // #pragma omp critical
  //       { ok = false; }
  //     }
  //   }
  //   assert(ok);
  //   std::cout << ok << std::endl;
  //   std::cout << "<<\n";
  // for (auto i : sa) std::cout << i << ' ';
  // std::cout << std::endl;
  // for (auto i : rank) std::cout << i << ' ';
  // std::cout << std::endl;
  // for (auto i : is_head) std::cout << (int)i << ' ';
  // std::cout << std::endl;
  // for (int i = 0; i < 50; i++) std::cout << rank[i] << ' ';
  // std::cout << std::endl;
  // for (int i = 0; i < 50; i++) std::cout << sa[i] << ' ';
  // std::cout << std::endl;
  sw2 = spdlog::stopwatch{};
  merge_sa_blocks(sa, rank, buf, is_head, index_offset);
  SPDLOG_DEBUG("merge sa blocks elapsed {}", sw2);
  // #pragma omp parallel for num_threads(NUM_THREADS)
  //   for (auto i = size_type{}; i < sa.size() - 1; i++) {
  //     if (sa[i] == sa[i + 1]) {
  // #pragma omp critical
  //       { ok = false; }
  //     }
  //   }
  //   assert(ok);
  //   std::cout << ok << std::endl;
  // std::cout << "<<\n";
  // for (auto i : sa) std::cout << i << ' ';
  // std::cout << std::endl;
  // for (auto i : rank) std::cout << i << ' ';
  // std::cout << std::endl;
  // for (auto i : is_head) std::cout << (int)i << ' ';
  // std::cout << std::endl;
  // for (int i = 0; i < 50; i++) std::cout << rank[i] << ' ';
  // std::cout << std::endl;
  // for (int i = 0; i < 50; i++) std::cout << sa[i] << ' ';
  // std::cout << std::endl;
  sw2 = spdlog::stopwatch{};
  calculate_new_rank_head(sa, rank, buf, is_head, index_offset);
  SPDLOG_DEBUG("calculate new head elapsed {}", sw2);
  // #pragma omp parallel for num_threads(NUM_THREADS)
  //   for (auto i = size_type{}; i < sa.size() - 1; i++) {
  //     if (sa[i] == sa[i + 1]) {
  // #pragma omp critical
  //       { ok = false; }
  //     }
  //   }
  //   assert(ok);
  //   std::cout << ok << std::endl;
  // std::cout << "<<\n";
  // for (auto i : sa) std::cout << i << ' ';
  // std::cout << std::endl;
  // for (auto i : rank) std::cout << i << ' ';
  // std::cout << std::endl;
  // for (auto i : is_head) std::cout << (int)i << ' ';
  // std::cout << std::endl;
  // for (int i = 0; i < 50; i++) std::cout << rank[i] << ' ';
  // std::cout << std::endl;
  // for (int i = 0; i < 50; i++) std::cout << sa[i] << ' ';
  // std::cout << std::endl;
  sw2 = spdlog::stopwatch{};
  compact(sa, rank, buf, is_head, starting_position, index_offset, sort_len);
  SPDLOG_DEBUG("compact elapsed {}", sw2);
  // #pragma omp parallel for num_threads(NUM_THREADS)
  //   for (auto i = size_type{}; i < sa.size() - 1; i++) {
  //     if (sa[i] == sa[i + 1]) {
  // #pragma omp critical
  //       { ok = false; }
  //     }
  //   }
  // std::cout << "<<\n";
  // for (auto i : sa) std::cout << i << ' ';
  // std::cout << std::endl;
  // for (auto i : rank) std::cout << i << ' ';
  // std::cout << std::endl;
  // for (auto i : is_head) std::cout << (int)i << ' ';
  // std::cout << std::endl;
  // for (int i = 0; i < 50; i++) std::cout << rank[i] << ' ';
  // std::cout << std::endl;
  // for (int i = 0; i < 50; i++) std::cout << sa[i] << ' ';
  // std::cout << std::endl;
  // assert(ok);
  // std::cout << ok << std::endl;
}

template<typename size_type>
void
get_overall_rank(const std::ranges::random_access_range auto& sa,
                 const std::ranges::random_access_range auto& rank,
                 const std::ranges::random_access_range auto& is_head) {
  std::vector<size_type> have_head(NUM_THREADS), head_tag(NUM_THREADS),
    head_sa_index(NUM_THREADS);
#pragma omp parallel for num_threads(NUM_THREADS)
  for (auto tid = size_type{}; tid < NUM_THREADS; tid++) {
    auto n_ = (size_type)sa.size() - 1;
    auto [L, R] = get_type_block_range(n_ + 1, NUM_THREADS, tid);
    bool local_have_head = false;
    size_type local_head_tag = size_type{}, local_head_sa_index = size_type{};
    for (auto i = L; i < R; i++) {
      if (is_head[i]) {
        local_have_head = true;
        local_head_tag = rank[sa[i]];
        local_head_sa_index = i;
      }
    }
#pragma omp critical
    {
      have_head[tid] = local_have_head;
      head_tag[tid] = local_head_tag;
      head_sa_index[tid] = local_head_sa_index;
    }
  }
  for (auto i = size_type{1}; i < NUM_THREADS; i++) {
    if (!have_head[i]) {
      head_tag[i] = head_tag[i - 1];
      head_sa_index[i] = head_sa_index[i - 1];
    }
  }
#pragma omp parallel for num_threads(NUM_THREADS)
  for (auto tid = size_type{}; tid < NUM_THREADS; tid++) {
    auto n_ = (size_type)sa.size() - 1;
    auto [L, R] = get_type_block_range(n_ + 1, NUM_THREADS, tid);
    size_type local_head_tag = (tid ? head_tag[tid - 1] : size_type{});
    size_type local_head_sa_index
      = (tid ? head_sa_index[tid - 1] : size_type{});
    for (auto i = L; i < R; i++) {
      if (is_head[i]) {
        local_head_tag = rank[sa[i]];
        local_head_sa_index = i;
      } else {
        rank[sa[i]] = local_head_tag + (i - local_head_sa_index);
      }
    }
  }
}

template<typename size_type>
void
get_overall_sa(std::ranges::random_access_range auto& sa,
               const std::ranges::random_access_range auto& rank) {
  auto n2 = (size_type)sa.size() - 1;
#pragma omp parallel for num_threads(NUM_THREADS)
  for (auto i = size_type{}; i < n2 + 1; i++) { sa[i] = -1; }
#pragma omp parallel for num_threads(NUM_THREADS)
  for (auto i = size_type{}; i < n2 + 1; i++) { sa[rank[i]] = i; }
}

// the prefix doubling step of the new kISS algorithm.
// rank = first idx of the same-value block
template<typename size_type>
void
prefix_doubling(std::ranges::random_access_range auto& sa,
                std::ranges::random_access_range auto& rank,
                std::ranges::random_access_range auto& buf,
                const std::ranges::random_access_range auto& starting_position,
                size_type sort_len) {
  // init rank and sa
  auto sw1 = spdlog::stopwatch{};
  auto n2 = (size_type)sa.size() - 1;
  init_sa<size_type>(sa, rank, buf);
  auto is_head = TypeVector(n2 + 1, BLOCK_ELEM_TYPE::NONHEAD);
  init_rank<size_type>(sa, rank, buf, is_head);
  SPDLOG_DEBUG("preparing elapsed {}", sw1);
  // std::cout << "<\n";
  // for (auto i : sa) std::cout << i << ' ';
  // std::cout << std::endl;
  // for (auto i : rank) std::cout << i << ' ';
  // std::cout << std::endl;
  // for (auto i : is_head) std::cout << (int)i << ' ';
  // std::cout << std::endl;
  // prefix doubling
  auto sa_dup = sa;
  for (uint64_t i = 1; i < sort_len; i *= 2) {
    radix_sort<size_type>(sa_dup, rank, buf, is_head, starting_position, i,
                          sort_len);
    // std::cout << "<\n";
    // for (auto i : sa_dup) std::cout << i << ' ';
    // std::cout << std::endl;
    // for (auto i : rank) std::cout << i << ' ';
    // std::cout << std::endl;
    // for (auto i : is_head) std::cout << (int)i << ' ';
    // std::cout << std::endl;
    // for (int i = 0; i < 50; i++) std::cout << rank[i] << ' ';
    // std::cout << std::endl;
    // for (int i = 0; i < 50; i++) std::cout << sa[i] << ' ';
    // std::cout << std::endl;
    if (sa_dup.size() <= 1)
      break;
  }
  // epilogue
  // for (int i = 0; i < 200; i++) std::cout << rank[i] << ' ';
  // std::cout << std::endl;
  // for (int i = 0; i < 200; i++) std::cout << sa[i] << ' ';
  // std::cout << std::endl;
  sa = std::ranges::subrange(std::begin(sa_dup), std::begin(sa_dup) + (n2 + 1));
  auto sw3 = spdlog::stopwatch{};
  get_overall_rank<size_type>(sa_dup, rank, is_head);
  // for (auto i : rank) std::cout << i << ' ';
  // std::cout << std::endl;
  // #pragma omp parallel for num_threads(NUM_THREADS)
  //   for (auto i = size_type{}; i < sa.size() - 1; i++) {
  //     if (sa[i] == sa[i + 1])
  //       assert(0);
  //   }
  get_overall_sa<size_type>(sa, rank);
  // for (int i : sa) std::cout << i << ' ';
  // std::cout << std::endl;
  // #pragma omp parallel for num_threads(NUM_THREADS)
  //   for (auto i = size_type{}; i < sa.size() - 1; i++) {
  //     if (sa[i] == sa[i + 1])
  //       assert(0);
  //   }
  SPDLOG_DEBUG("epilogue elapsed {}", sw3);
}

template<typename size_type>
void
place_back_lms(const std::ranges::random_access_range auto& T,
               const std::ranges::random_access_range auto& sa,
               std::ranges::random_access_range auto& sorted_lms,
               const std::ranges::random_access_range auto& buf,
               const std::ranges::random_access_range auto& starting_position,
               const std::ranges::random_access_range auto& valid_position) {
  //   // pre scan
  //   auto n = (size_type)T.size();
  auto n2 = (size_type)sa.size() - 1;
  //   std::vector<size_type> have_lms(NUM_THREADS), first_lms(NUM_THREADS),
  //     last_lms(NUM_THREADS), len(NUM_THREADS);
  // #pragma omp parallel for num_threads(NUM_THREADS)
  //   for (auto tid = size_type{}; tid < NUM_THREADS; tid++) {
  //     auto [L, R] = get_type_block_range(n + 1, NUM_THREADS, tid);
  //     bool local_have_lms = false;
  //     size_type local_first_lms = 0, local_last_lms = 0;
  //     size_type local_len = 0;
  //     size_type last_lms_index = 0;
  //     for (auto i = L; i < R; i++) {
  //       if (is_LMS(T, i)) {
  //         if (!local_have_lms) {
  //           local_first_lms = local_last_lms = i;
  //         } else {
  //           local_last_lms = i;
  //           local_len += (i - last_lms_index + compress_block_length)
  //                        / compress_block_length;
  //         }
  //         local_have_lms = true;
  //         last_lms_index = i;
  //       }
  //     }
  // #pragma omp critical
  //     {
  //       have_lms[tid] = local_have_lms;
  //       first_lms[tid] = local_first_lms;
  //       last_lms[tid] = local_last_lms;
  //       len[tid] = local_len;
  //     }
  //   }
  //   // prefix sum to get number of char before each block
  //   bool have_last_block_lms = false;
  //   size_type last_block_lms = size_type{};
  //   for (auto i = size_type{}; i < NUM_THREADS; i++) {
  //     if (have_lms[i]) {
  //       if (have_last_block_lms) {
  //         len[i - 1] += (first_lms[i] - last_block_lms +
  //         compress_block_length)
  //                       / compress_block_length;
  //       }
  //       have_last_block_lms = true;
  //       last_block_lms = last_lms[i];
  //     }
  //   }
  //   std::inclusive_scan(std::begin(len), std::end(len), std::begin(len));
  //   // do mapping: n2 index -> normal index (+ boolean array to know whether
  //   each
  //   // index is valid)
  //   TypeVector valid(n2, 0);
  // #pragma omp parallel for num_threads(NUM_THREADS)
  //   for (auto tid = size_type{}; tid < NUM_THREADS; tid++) {
  //     auto [L, R] = get_type_block_range(n + 1, NUM_THREADS, tid);
  //     auto index = (tid ? len[tid - 1] : size_type{});
  //     bool local_have_lms = false;
  //     auto last_lms_index = L;
  //     for (auto i = L; i < R; i++) {
  //       if (is_LMS(T, i)) {
  //         if (local_have_lms) {
  //           index += (i - last_lms_index + compress_block_length)
  //                    / compress_block_length;
  //         }
  //         valid[index] = true;
  //         buf[index] = i;
  //         local_have_lms = true;
  //         last_lms_index = i;
  //       }
  //     }
  //   }
  // scan sa to put them in correct order
  std::vector<size_type> count_sorted_lms(NUM_THREADS);
#pragma omp parallel for num_threads(NUM_THREADS)
  for (auto tid = size_type{}; tid < NUM_THREADS; tid++) {
    auto [L, R] = get_type_block_range(n2 + 1, NUM_THREADS, tid);
    auto local_count_sorted_lms = size_type{};
    for (auto i = L; i < R; i++) {
      if (valid_position[sa[i]])
        local_count_sorted_lms++;
    }
#pragma omp critical
    { count_sorted_lms[tid] = local_count_sorted_lms; }
  }
  std::inclusive_scan(std::begin(count_sorted_lms), std::end(count_sorted_lms),
                      std::begin(count_sorted_lms));
#pragma omp parallel for num_threads(NUM_THREADS)
  for (auto tid = size_type{}; tid < NUM_THREADS; tid++) {
    auto [L, R] = get_type_block_range(n2 + 1, NUM_THREADS, tid);
    auto local_count_sorted_lms
      = (tid ? count_sorted_lms[tid - 1] : size_type{});
    for (auto i = L; i < R; i++) {
      if (valid_position[sa[i]])
        sorted_lms[local_count_sorted_lms++] = starting_position[sa[i]];
    }
  }
}
#undef MASK
#undef NUM_THREADS
#undef PREFIX_DOUBLING_NUM_THREADS
}  // namespace kiss

}  // namespace biovoltron