// FIXME: refactor the code!
#pragma once

#include <spdlog/spdlog.h>
#include <spdlog/stopwatch.h>

#include <algorithm>
#include <biovoltron/container/xbit_vector.hpp>
#include <biovoltron/utility/istring.hpp>
#include <biovoltron/utility/thread_pool.hpp>
#include <bit>
#include <chrono>
#include <cstdlib>
#include <execution>
#include <iostream>
#include <numeric>
#include <tuple>
#include <vector>
#include <xmmintrin.h>

#include "omp.h"

namespace biovoltron {

namespace kiss {

#define MASK        ((1 << 16) - 1)
#define NUM_THREADS std::thread::hardware_concurrency()
#define PREFIX_DOUBLING_NUM_THREADS(n) \
  std::min((uint32_t)NUM_THREADS, (uint32_t)n)
#define LARGE_SEGMENT_THRESHOLD (1 << 14)

enum SUFFIX_TYPE { L_TYPE = 0, S_TYPE = 1 };
enum BLOCK_ELEM_TYPE { NONHEAD = 0, HEAD = 1 };
using TypeVector = biovoltron::detail::XbitVector<1, std::uint8_t,
                                                  std::allocator<std::uint8_t>>;

// is T[i] an LMS character?
template<typename size_type>
auto
is_LMS(const auto& T, size_type i) {
  auto n = (size_type)T.size();
  return i == n
         or (i > 0 and T[i - 1] == SUFFIX_TYPE::L_TYPE
             and T[i] == SUFFIX_TYPE::S_TYPE);
}

// Partition num_items of items into num_blocks blocks, the function
// outputs the range for block number block_idx.
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

// Get the length of the encoded string given the compress_block_length (l).
template<typename size_type>
auto
get_encoded_reference_length(const auto& lms, size_type l) {
  auto n1 = (size_type)lms.size();
  auto length = size_type{};
#pragma omp parallel for num_threads(NUM_THREADS)
  for (auto tid = size_type{}; tid < NUM_THREADS; tid++) {
    auto [L, R] = get_type_block_range(n1, NUM_THREADS, tid);
    auto local_length_count = size_type{};
    for (auto i = L; i + 1 < R; i++) {
      local_length_count += (lms[i + 1] - lms[i] + l - 1) / l;
    }
    if (L != R && R < n1)
      local_length_count += (lms[R] - lms[R - 1] + l - 1) / l;
#pragma omp atomic
    length += local_length_count;
  }
  return length;
}

// Reduce the string to encoded form by assigning every
// l-mer into an integer.
template<typename size_type>
void
encode_reference(const std::ranges::random_access_range auto& S, auto& lms,
                 auto& cS, const std::ranges::random_access_range auto& buf,
                 const std::ranges::random_access_range auto& starting_position,
                 std::ranges::random_access_range auto& valid_position,
                 size_type K, size_type l) {
  auto n = (size_type)S.size();
  auto n1 = (size_type)lms.size();
  auto n2 = (size_type)cS.size();
  std::vector<size_type> block_cS_start_index(NUM_THREADS);
  // Find the starting index of cS for each thread
#pragma omp parallel for num_threads(NUM_THREADS)
  for (auto tid = size_type{}; tid < NUM_THREADS; tid++) {
    auto [L, R] = get_type_block_range(n1, NUM_THREADS, tid);
    auto local_block_cS_start_index = size_type{};
    for (auto i = L; i + 1 < R; i++) {
      if (i + 1 < n1)
        local_block_cS_start_index += (lms[i + 1] - lms[i] + l - 1) / l;
    }
    if (L != R && R < n1)
      local_block_cS_start_index += (lms[R] - lms[R - 1] + l - 1) / l;
#pragma omp critical
    block_cS_start_index[tid] = local_block_cS_start_index;
  }
  // Exclusive scan
  std::exclusive_scan(block_cS_start_index.begin(), block_cS_start_index.end(),
                      block_cS_start_index.begin(), size_type{}, std::plus<>{});
  // Fill in the encoded value
  auto bits_for_character = std::bit_width(K - 1);
#pragma omp parallel for num_threads(NUM_THREADS)
  for (auto tid = size_type{}; tid < NUM_THREADS; tid++) {
    auto [L, R] = get_type_block_range(n1, NUM_THREADS, tid);
    auto local_cS_index = block_cS_start_index[tid];
    for (auto i = L; i < R; i++) {
      if (i + 1 >= n1)
        continue;
      auto encoded_LMS_substring_length = (lms[i + 1] - lms[i] + l - 1) / l;
      auto l_border = lms[i];
      auto r_border = l_border + l;
      for (auto j = 0; j < encoded_LMS_substring_length; j++) {
        auto encoded = size_type{};
        for (auto k = l_border; k < r_border; k++) {
          auto character = (k < n ? S[k] : 0);
          encoded = (encoded << bits_for_character) | character;
        }
        cS[local_cS_index] = encoded;
        starting_position[local_cS_index] = l_border;
        valid_position[local_cS_index] = (j == 0);
        l_border += l;
        r_border += l;
        local_cS_index++;
      }
    }
  }
  starting_position[n2] = n;
}

template<typename size_type>
auto
get_key(const std::ranges::random_access_range auto& arr, size_type idx,
        size_type index_offset) {
  return ((uint64_t)idx + index_offset >= arr.size()) ?
           size_type{} :
           arr[idx + index_offset];
}

template<typename size_type>
void
sort_range(std::ranges::random_access_range auto& arr,
           std::ranges::random_access_range auto& buf,
           auto&& member, size_type bits, size_type index_offset,
           size_type num_threads) {
  if (bits & 0xF) {
    throw std::invalid_argument("Invalid bits argument");
    return;
  }
  if (arr.size() != buf.size()) {
    throw std::invalid_argument("arr and buf size mismatch");
    return;
  }

  auto n = (size_type)arr.size();

  if (n <= LARGE_SEGMENT_THRESHOLD) {
    size_type bitmask = (size_type{1} << bits) - 1;
    auto comparator = [&member, index_offset, bitmask](size_type i, size_type j) {
      return (get_key(member, i, index_offset) & bitmask) 
           < (get_key(member, j, index_offset) & bitmask);
    };
    std::sort(std::begin(arr), std::end(arr), comparator);
    return;
  }

  auto passes = bits / 16;

  for (auto k = size_type{}; k < passes; k++) {
      auto offsets = std::vector(num_threads, std::vector(1 << 16, size_type{}));
      auto right_shift_offset = k * 16;
    // calculate count for each mask of each thread
#pragma omp parallel for num_threads(num_threads)
    for (auto i = size_type{}; i < arr.size(); i++) {
      auto tid = omp_get_thread_num();
      auto& elem = ((k & 1) ? buf[i] : arr[i]);
      auto key = (member[elem] >> right_shift_offset) & (0xFFFF);
      offsets[tid][key]++;
    }
    // Calculate offset for each mask of each thread
    auto total = std::vector<size_type>(1 << 16, 0);
#pragma omp parallel for num_threads(num_threads)
    for (auto i = 0; i < (1 << 16); i++) {
      auto& x = total[i];
      for (auto tid = 0; tid < num_threads; tid++) {
        x += offsets[tid][i];
        offsets[tid][i] = x - offsets[tid][i];
      }
    }
    std::exclusive_scan(total.begin(), total.end(), total.begin(), size_type{},
                        std::plus<>{});
    // Place back the elements
  #pragma omp parallel for num_threads(num_threads)
    for (auto i = size_type{}; i < arr.size(); i++) {
      auto tid = omp_get_thread_num();
      auto& elem = ((k & 1) ? buf[i] : arr[i]);
      auto key = (member[elem] >> right_shift_offset) & (0xFFFF);
      auto idx = total[key] + (offsets[tid][key]++);
      auto& output_addr = ((k & 1) ? arr[idx] : buf[idx]);
      output_addr = elem;
    }
  }
  if (passes & 1) {
    std::memmove(arr.data(), buf.data(), (uint64_t)n * sizeof(size_type));
  }
}

// Perform 1-sort on sa.
template<typename size_type>
void
initialize_sa(const std::ranges::random_access_range auto& sa,
              const std::ranges::random_access_range auto& rank,
              const std::ranges::random_access_range auto& buf) {
  auto n2 = (size_type)sa.size();
#pragma omp parallel for num_threads(NUM_THREADS)
  for (auto i = size_type{}; i < n2; i++) { sa[i] = i; }
  std::swap(sa[n2 - 1], sa[0]);  // to make sure that the sentinel would be
                                 // the first
  // perform radix sort based on key rank[i]
  sort_range<size_type>(sa, buf, rank, sizeof(size_type) * 8, 0, NUM_THREADS);
}

// Calculate the inverse suffix array.
template<typename size_type>
void
initialize_rank(const std::ranges::random_access_range auto& sa,
                std::ranges::random_access_range auto& rank,
                std::ranges::random_access_range auto& buf,
                std::ranges::random_access_range auto& is_head) {
  auto n2 = (size_type)sa.size();
  std::vector<bool> is_block_has_head(NUM_THREADS);
  std::vector<size_type> block_last_head(NUM_THREADS);
  // pre scan
#pragma omp parallel for num_threads(NUM_THREADS)
  for (auto tid = size_type{}; tid < NUM_THREADS; tid++) {
    auto [L, R] = get_type_block_range(n2, NUM_THREADS, tid);
    auto local_block_last_head = L;
    auto local_block_accumulated_length = size_type{};
    bool local_block_has_head = false;
    for (auto i = L; i < R; i++) {
      size_type is_head_i = (i == 0 || sa[i - 1] == n2 - 1 || rank[sa[i]] != rank[sa[i - 1]]);
      is_head[i] = is_head_i;
      local_block_has_head = (local_block_has_head || is_head_i);
      local_block_accumulated_length += 1;
      local_block_last_head += local_block_accumulated_length & (-is_head_i);
      buf[sa[i]] = local_block_last_head;
      local_block_accumulated_length &= (~(-is_head_i));
    }
#pragma omp critical
    {
      is_block_has_head[tid] = local_block_has_head;
      block_last_head[tid] = local_block_last_head;
    }
  }
  // collect last labels
  for (auto tid = size_type{1}; tid < NUM_THREADS; tid++) {
    if (!is_block_has_head[tid]) {
      block_last_head[tid] = block_last_head[tid - 1];
    }
  }
  // post scan
#pragma omp parallel for num_threads(NUM_THREADS)
  for (auto tid = size_type{}; tid < NUM_THREADS; tid++) {
    auto [L, R] = get_type_block_range(n2, NUM_THREADS, tid);
    auto local_block_last_head = (tid > 0 ? block_last_head[tid - 1] : 0);
    auto i = L;
    while (i < R && !is_head[i]) {
      buf[sa[i]] = local_block_last_head;
      ++i;
    }
  }
  // extra swap
  auto new_rank = std::ranges::subrange(std::begin(buf), std::begin(buf) + n2);
  buf = rank;
  rank = new_rank;
}

// Radix sort on range [L, R] in sa by key.
template<typename size_type>
void
sort_sa_same_first_value(const std::ranges::random_access_range auto& sa,
                         const std::ranges::random_access_range auto& rank,
                         const std::ranges::random_access_range auto& buf,
                         size_type index_offset, size_type L, size_type R) {
  auto comparator = [&rank, index_offset](size_type i, size_type j) {
    return get_key(rank, i, index_offset) < get_key(rank, j, index_offset);
  };
  std::sort(std::begin(sa) + L, std::begin(sa) + R, comparator);
}

// Sort each segment in the same block.
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
        sort_sa_same_first_value(sa, rank, buf, index_offset, last_left, i);
        last_left = i;
      }
    }
    sort_sa_same_first_value(sa, rank, buf, index_offset, last_left, R);
  }
}

// TODO: more details
// The recorded data for each block in SegmentedSort.
// block_has_head: whether the block contains head of a segment
// first_segment_end: the ending index of the first segment
// last_segment_head: the starting index of the last segment
template<typename size_type>
struct MergeBlockData {
  bool is_block_has_head;
  size_type first_segment_end;
  size_type last_segment_head;
  MergeBlockData()
  : is_block_has_head(false),
    first_segment_end(size_type{}),
    last_segment_head(size_type{}) { }
  MergeBlockData(bool is_block_has_head_, size_type first_segment_end_,
                 size_type last_segment_head_)
  : is_block_has_head(is_block_has_head_),
    first_segment_end(first_segment_end_),
    last_segment_head(last_segment_head_) { }
};

template<typename size_type>
MergeBlockData<size_type>
merge_sa_blocks_recursive(
  const std::ranges::random_access_range auto& sa,
  const std::ranges::random_access_range auto& rank,
  const std::ranges::random_access_range auto& buf,
  std::vector<MergeBlockData<size_type>>& merge_block_data_list,
  size_type index_offset, size_type L, size_type R) {
  if (L + 1 == R) {
    return merge_block_data_list[L];
  }
  auto mid = (L + R) / 2;
  auto [l_has_head, l_first_segment_end, l_last_segment_head]
    = merge_sa_blocks_recursive(sa, rank, buf, merge_block_data_list,
                                index_offset, L, mid);
  auto [r_has_head, r_first_segment_end, r_last_segment_head]
    = merge_sa_blocks_recursive(sa, rank, buf, merge_block_data_list,
                                index_offset, mid, R);
  auto n_ = (size_type)sa.size() - 1;
  auto [mid_block_start, _] = get_type_block_range(n_ + 1, NUM_THREADS, mid);
  bool is_indicies_in_range
    = (l_last_segment_head < n_ + 1 && mid_block_start < n_ + 1
       && r_first_segment_end <= n_ + 1);
  bool is_ranges_in_order = (l_last_segment_head < mid_block_start
                             && mid_block_start < r_first_segment_end);
  if (is_indicies_in_range && is_ranges_in_order) {
    auto comparator = [&](size_type idx1, size_type idx2) {
      return get_key(rank, idx1, index_offset)
             < get_key(rank, idx2, index_offset);
    };
    std::merge(std::begin(sa) + l_last_segment_head,
               std::begin(sa) + mid_block_start,
               std::begin(sa) + mid_block_start,
               std::begin(sa) + r_first_segment_end,
               std::begin(buf) + l_last_segment_head, comparator);
    std::swap_ranges(std::begin(sa) + l_last_segment_head,
                     std::begin(sa) + r_first_segment_end,
                     std::begin(buf) + l_last_segment_head);
  }
  auto new_has_head = (l_has_head || r_has_head);
  auto new_first_segment_end
    = (l_has_head ? l_first_segment_end : r_first_segment_end);
  auto new_last_segment_head
    = (r_has_head ? r_last_segment_head : l_last_segment_head);
  return MergeBlockData(new_has_head, new_first_segment_end,
                        new_last_segment_head);
}

// merge the blocks by merging the consecutive segments.
template<typename size_type>
void
merge_sa_blocks(const std::ranges::random_access_range auto& sa,
                const std::ranges::random_access_range auto& rank,
                const std::ranges::random_access_range auto& buf,
                const std::ranges::random_access_range auto& is_head,
                size_type index_offset) {
  std::vector<MergeBlockData<size_type>> merge_block_data_list(NUM_THREADS);
#pragma omp parallel for num_threads(NUM_THREADS)
  for (auto tid = size_type{}; tid < NUM_THREADS; tid++) {
    auto n_ = (size_type)sa.size() - 1;
    auto [L, R] = get_type_block_range(n_ + 1, NUM_THREADS, tid);

    size_type local_first_segment_end = L, local_last_segment_head = L;
    bool have_true = false;
    for (size_type i = L; i < R; i++) {
      if (!have_true)
        local_first_segment_end = i;
      if (is_head[i]) {
        have_true = true;
        local_last_segment_head = i;
      }
    }
    if (!have_true && L < n_ + 1)
      local_first_segment_end++;
#pragma omp critical
    {
      merge_block_data_list[tid] = MergeBlockData(
        have_true, local_first_segment_end, local_last_segment_head);
    }
  }
  merge_sa_blocks_recursive(sa, rank, buf, merge_block_data_list, index_offset,
                            size_type{}, NUM_THREADS);
}

template<typename size_type>
void
calculate_new_rank_head(const std::ranges::random_access_range auto& sa,
                        std::ranges::random_access_range auto& rank,
                        std::ranges::random_access_range auto& is_head,
                        size_type index_offset) {
  auto n_ = (size_type)sa.size();
  auto is_new_head = TypeVector(n_, SUFFIX_TYPE::L_TYPE);
  std::vector<size_type> head_tag(NUM_THREADS), head_sa_index(NUM_THREADS);
  std::vector<int> level(NUM_THREADS);
#pragma omp parallel for num_threads(NUM_THREADS)
  for (auto tid = size_type{}; tid < NUM_THREADS; tid++) {
    auto [L, R] = get_type_block_range(n_, NUM_THREADS, tid);
    size_type local_head_tag = size_type{}, local_head_sa_index = size_type{};
    int local_level = 0;
    for (auto i = L; i < R; i++) {
      if (is_head[i]) {
        local_level = 2;
        local_head_tag = rank[sa[i]];
        local_head_sa_index = i;
      } else if (get_key<size_type>(rank, sa[i], index_offset)
                 != get_key<size_type>(rank, sa[i - 1], index_offset)) {
        is_new_head[i] = 1;
        local_level = std::max(local_level, 1);
        local_head_tag = (i - local_head_sa_index) + local_head_tag;
        local_head_sa_index = i;
      }
    }
#pragma omp critical
    {
      head_tag[tid] = local_head_tag;
      head_sa_index[tid] = local_head_sa_index;
      level[tid] = local_level;
    }
  }
  for (auto i = size_type{1}; i < NUM_THREADS; i++) {
    if (level[i] == 0) {
      head_tag[i] = head_tag[i - 1];
      head_sa_index[i] = head_sa_index[i - 1];
    } else if (level[i] == 1) {
      head_tag[i] = head_tag[i - 1] + (head_sa_index[i] - head_sa_index[i - 1]);
    }
  }
#pragma omp parallel for num_threads(NUM_THREADS)
  for (auto tid = size_type{}; tid < NUM_THREADS; tid++) {
    auto [L, R] = get_type_block_range(n_, NUM_THREADS, tid);
    size_type local_head_tag = (tid ? head_tag[tid - 1] : size_type{});
    size_type local_head_sa_index
      = (tid ? head_sa_index[tid - 1] : size_type{});
    for (auto i = L; i < R; i++) {
      if (is_head[i]) {
        local_head_tag = rank[sa[i]];
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
}

template<typename size_type>
void
compact(std::ranges::random_access_range auto& sa,
        std::ranges::random_access_range auto& buf,
        std::ranges::random_access_range auto& is_head) {
  auto n_ = (size_type)sa.size();
  auto n2 = (size_type)buf.size();
  std::vector<size_type> count_compact(NUM_THREADS + 1);
#pragma omp parallel for num_threads(NUM_THREADS)
  for (auto tid = size_type{}; tid < NUM_THREADS; tid++) {
    size_type local_count_compact = size_type{};
    auto [L, R] = get_type_block_range(n_, NUM_THREADS, tid);
    for (auto i = L; i < R; i++) {
      bool singleton = is_head[i] && (i + 1 >= n_ || is_head[i + 1]);
      if (!singleton)
        local_count_compact++;
    }
#pragma omp critical
    { count_compact[tid + 1] = local_count_compact; }
  }
  std::inclusive_scan(std::begin(count_compact), std::end(count_compact),
                      std::begin(count_compact));
  auto new_n_ = count_compact[NUM_THREADS];
  auto is_head_buf = TypeVector(new_n_, BLOCK_ELEM_TYPE::NONHEAD);
#pragma omp parallel for num_threads(NUM_THREADS)
  for (auto tid = size_type{}; tid < NUM_THREADS; tid++) {
    size_type local_count_compact = count_compact[tid];
    auto [L, R] = get_type_block_range(n_, NUM_THREADS, tid);
    auto newL = count_compact[tid], newR = count_compact[tid + 1];
    for (auto i = L; i < R; i++) {
      bool singleton = is_head[i] && (i + 1 >= n_ || is_head[i + 1]);
      if (!singleton) {
        buf[local_count_compact] = sa[i];
        if (local_count_compact - newL < 8
            || newR - local_count_compact
                 < 8) {  // prevent race condition in is_head_buf
#pragma omp critical
          { is_head_buf[local_count_compact] = is_head[i]; }
        } else {
          is_head_buf[local_count_compact] = is_head[i];
        }
        local_count_compact++;
      }
    }
  }
  auto new_buf = std::ranges::subrange(std::begin(sa), std::begin(sa) + (n2));
  sa = buf;
  buf = new_buf;
  sa = std::ranges::subrange(std::begin(sa), std::begin(sa) + (new_n_));
  is_head = is_head_buf;
}

template<typename size_type>
void
h_sort(std::ranges::random_access_range auto& sa,
       std::ranges::random_access_range auto& rank,
       std::ranges::random_access_range auto& buf,
       std::ranges::random_access_range auto& is_head, size_type index_offset) {
  auto sw2 = spdlog::stopwatch{};
  sort_sa_blocks(sa, rank, buf, is_head, index_offset);
  SPDLOG_DEBUG("sort sa blocks elapsed {}", sw2);
  sw2 = spdlog::stopwatch{};
  merge_sa_blocks(sa, rank, buf, is_head, index_offset);
  SPDLOG_DEBUG("merge sa blocks elapsed {}", sw2);
  sw2 = spdlog::stopwatch{};
  calculate_new_rank_head(sa, rank, is_head, index_offset);
  SPDLOG_DEBUG("calculate new head elapsed {}", sw2);
  sw2 = spdlog::stopwatch{};
  compact<size_type>(sa, buf, is_head);
  SPDLOG_DEBUG("compact elapsed {}", sw2);
}

template<typename size_type>
void
get_overall_rank(const std::ranges::random_access_range auto& sa,
                 const std::ranges::random_access_range auto& rank,
                 const std::ranges::random_access_range auto& is_head) {
  std::vector<size_type> block_has_head(NUM_THREADS), head_tag(NUM_THREADS),
    head_sa_index(NUM_THREADS);
  auto n_ = (size_type)sa.size();
#pragma omp parallel for num_threads(NUM_THREADS)
  for (auto tid = size_type{}; tid < NUM_THREADS; tid++) {
    auto [L, R] = get_type_block_range(n_, NUM_THREADS, tid);
    bool local_block_has_head = false;
    size_type local_head_tag = size_type{}, local_head_sa_index = size_type{};
    for (auto i = L; i < R; i++) {
      if (is_head[i]) {
        local_block_has_head = true;
        local_head_tag = rank[sa[i]];
        local_head_sa_index = i;
      }
    }
#pragma omp critical
    {
      block_has_head[tid] = local_block_has_head;
      head_tag[tid] = local_head_tag;
      head_sa_index[tid] = local_head_sa_index;
    }
  }
  for (auto i = size_type{1}; i < NUM_THREADS; i++) {
    if (!block_has_head[i]) {
      head_tag[i] = head_tag[i - 1];
      head_sa_index[i] = head_sa_index[i - 1];
    }
  }
#pragma omp parallel for num_threads(NUM_THREADS)
  for (auto tid = size_type{}; tid < NUM_THREADS; tid++) {
    auto [L, R] = get_type_block_range(n_, NUM_THREADS, tid);
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
  auto n2 = (size_type)sa.size();
#pragma omp parallel for num_threads(NUM_THREADS)
  for (auto i = size_type{}; i < n2; i++) { sa[rank[i] - 1] = i; }
}

// the prefix doubling step of the new kISS algorithm.
// rank = first idx of the same-value block
template<typename size_type>
void
prefix_doubling(std::ranges::random_access_range auto& sa,
                std::ranges::random_access_range auto& rank,
                std::ranges::random_access_range auto& buf,
                size_type sort_len) {
  // init rank and sa
  auto sw1 = spdlog::stopwatch{};
  auto n2 = (size_type)sa.size();
  initialize_sa<size_type>(sa, rank, buf);
  SPDLOG_DEBUG("preparing 1 elapsed {}", sw1);
  sw1 = spdlog::stopwatch{};
  // This Boolean array is used for easier implementation.
  auto is_head = TypeVector(n2, BLOCK_ELEM_TYPE::NONHEAD);
  initialize_rank<size_type>(sa, rank, buf, is_head);
  SPDLOG_DEBUG("preparing 2 elapsed {}", sw1);
  auto sa_dup = sa;
  for (uint64_t i = 1; i < sort_len; i *= 2) {
    h_sort<size_type>(sa_dup, rank, buf, is_head, i);
    if (sa_dup.size() <= 1)
      break;
  }
  // epilogue
  auto sw3 = spdlog::stopwatch{};
  get_overall_rank<size_type>(sa_dup, rank, is_head);
  sa = std::ranges::subrange(std::begin(sa_dup),
                             std::begin(sa_dup) + (sa.size()));
  get_overall_sa<size_type>(sa, rank);
  SPDLOG_DEBUG("epilogue elapsed {}", sw3);
}

template<typename size_type>
void
place_back_lms(const std::ranges::random_access_range auto& sa,
               std::ranges::random_access_range auto& sorted_lms,
               const std::ranges::random_access_range auto& starting_position,
               const std::ranges::random_access_range auto& valid_position) {
  auto n2 = (size_type)sa.size();
  // scan sa to put them in correct order
  std::vector<size_type> count_sorted_lms(NUM_THREADS);
#pragma omp parallel for num_threads(NUM_THREADS)
  for (auto tid = size_type{}; tid < NUM_THREADS; tid++) {
    auto [L, R] = get_type_block_range(n2, NUM_THREADS, tid);
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
    auto [L, R] = get_type_block_range(n2, NUM_THREADS, tid);
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
