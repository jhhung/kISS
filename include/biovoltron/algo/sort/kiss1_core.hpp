#pragma once

#include <biovoltron/algo/sort/kiss_common.hpp>
#include <biovoltron/algo/sort/structs.hpp>
#include <biovoltron/algo/sort/utils.hpp>
#include <spdlog/spdlog.h>
#include <spdlog/stopwatch.h>

#include <omp.h>
#include <thread>
#include <ranges>

namespace biovoltron {

namespace kiss {

template <typename char_type, typename size_type>
void lms_suffix_direct_sort_dna(
  const std::ranges::random_access_range auto &S,
  const std::ranges::random_access_range auto &SA,
  const std::ranges::random_access_range auto &buffer,
  size_type m,
  size_type k,
  const size_t num_threads,
  vector<ThreadState<size_type>> &states
) {
  auto packed_S = PackedDNAString<char_type, size_type>(S, states, num_threads);

#pragma omp parallel num_threads(num_threads)
  {
    auto tid = omp_get_thread_num();
    auto &state = states[tid];
    distribute_workload<size_type>(m, state);

    for (auto i = state.beg; i < state.beg + state.len; i++) {
      auto prefix = packed_S.load_prefix_length_less_than_16(SA[i], KISS1_SPLIT_SORT_PREFIX_SIZE);
      state.buffer[prefix]++;
    }
  }

  std::array<size_type, KISS1_SPLIT_SORT_BUCKET_SIZE + 1> segment_start{};

#pragma omp parallel for num_threads(num_threads)
  for (auto i = size_type{}; i < KISS1_SPLIT_SORT_BUCKET_SIZE; i++) {
    for (auto j = size_t{}; j < num_threads; j++) {
      segment_start[i] += states[j].buffer[i];
      states[j].buffer[i] = segment_start[i] - states[j].buffer[i];
    }
  }

  std::exclusive_scan(std::begin(segment_start), std::end(segment_start),
      std::begin(segment_start), size_type{});

#pragma omp parallel num_threads(num_threads)
  {
    auto tid = omp_get_thread_num();
    auto &state = states[tid];
    distribute_workload<size_type>(m, state);

    for (auto i = state.beg; i < state.beg + state.len; i++) {
      auto prefix = packed_S.load_prefix_length_less_than_16(SA[i], KISS1_SPLIT_SORT_PREFIX_SIZE);
      auto idx = segment_start[prefix] + (state.buffer[prefix]++);
      buffer[idx] = SA[i];
    }
  }

  auto n = S.size();

  auto comparator = [&S, k, n, &packed_S](size_type i, size_type j) {
    auto sorted_len = size_type{};
    while (sorted_len <= k && 
           i + KISS1_SPLIT_SORT_STRIDE_DNA <= n &&
           j + KISS1_SPLIT_SORT_STRIDE_DNA <= n) {
      if (i + KISS1_SPLIT_SORT_STRIDE_DNA > n)
        break;
      if (j + KISS1_SPLIT_SORT_STRIDE_DNA > n)
        break; 

      auto a = packed_S.load_prefix_length_125(i);
      auto b = packed_S.load_prefix_length_125(j);
      auto eq = _mm256_cmpeq_epi8(a, b);
      auto neq_mask = ~((uint32_t)_mm256_movemask_epi8(eq));

      if (neq_mask != 0) {
        auto msb_mask = (1UL << (31 - _lzcnt_u32(neq_mask)));
        auto gt = _mm256_max_epu8(a, b);
        gt = _mm256_cmpeq_epi8(gt, a);
        auto gt_mask = (uint32_t)_mm256_movemask_epi8(gt);
        return !((msb_mask & gt_mask) > 0);
      }
      
      sorted_len += KISS1_SPLIT_SORT_STRIDE_DNA;
      i += KISS1_SPLIT_SORT_STRIDE_DNA;
      j += KISS1_SPLIT_SORT_STRIDE_DNA;
    }
    while (sorted_len + 1 <= k && 
           i < n && 
           j < n) {
      if (S[i] != S[j])
        return S[i] < S[j];

      sorted_len++;
      i++;
      j++;
    }
    if (sorted_len >= k) {
      return i < j;
    }
    return (i == n);
  };

#pragma omp parallel for schedule(dynamic) num_threads(num_threads) 
  for (auto i = size_type{}; i < KISS1_SPLIT_SORT_BUCKET_SIZE; i++) {
    std::sort(std::begin(buffer) + segment_start[i], 
              std::begin(buffer) + segment_start[i + 1],
              comparator);
  }

  std::memmove(&SA[0], &buffer[0], sizeof(size_type) * m);
}

template <typename char_type, typename size_type>
void lms_suffix_direct_sort(
  const std::ranges::random_access_range auto &S,
  vector<size_type> &SA,
  size_type m,
  size_type k
) {
  auto n = S.size();

  auto comparator = [&S, k, n](size_type i, size_type j) { // {{{
    auto sorted_len = size_type{};
    while (sorted_len <= k && 
           i + KISS1_SPLIT_SORT_STRIDE <= n &&
           j + KISS1_SPLIT_SORT_STRIDE <= n) {
      if (i + KISS1_SPLIT_SORT_STRIDE > n)
        break;
      if (j + KISS1_SPLIT_SORT_STRIDE > n)
        break; 

      auto a = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(&S[i]));
      auto b = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(&S[j]));
      auto eq = _mm256_cmpeq_epi8(a, b);
      auto neq_mask = ~((uint32_t)_mm256_movemask_epi8(eq));

      if (neq_mask != 0) {
        auto lsb_mask = (1UL << _tzcnt_u32(neq_mask));
        auto gt = _mm256_max_epu8(a, b);
        gt = _mm256_cmpeq_epi8(gt, a);
        auto gt_mask = (uint32_t)_mm256_movemask_epi8(gt);
        return !((lsb_mask & gt_mask) > 0);
      }

      sorted_len += KISS1_SPLIT_SORT_STRIDE;
      i += KISS1_SPLIT_SORT_STRIDE;
      j += KISS1_SPLIT_SORT_STRIDE;
    }
    while (sorted_len + 1 <= k && 
           i < n && 
           j < n) {
      if (S[i] != S[j])
        return S[i] < S[j];

      sorted_len++;
      i++;
      j++;
    }
    if (sorted_len >= k) {
      return i < j;
    }
    return (i == n);
  }; // }}}

  auto small_cmp = [&S, k](auto i, auto j) {
    auto ri = std::ranges::subrange(S.begin() + i, S.begin() + std::min<size_type>(S.size(), i + k));
    auto rj = std::ranges::subrange(S.begin() + j, S.begin() + std::min<size_type>(S.size(), j + k));
    return std::ranges::lexicographical_compare(ri, rj);
  };

  const size_type thres = 100'000'000;
  if (n < thres) {
    std::sort(std::execution::par_unseq,
              std::begin(SA),
              std::begin(SA) + m,
              small_cmp);
  } else {
    std::sort(std::execution::par_unseq,
              std::begin(SA),
              std::begin(SA) + m,
              comparator);
  }
}


template <typename char_type, typename size_type>
void kiss1_suffix_array_dna(
  const vector<char_type> &S,
  vector<size_type> &SA,
  size_type k = 256u,
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
  auto buffer = std::ranges::subrange(std::begin(SA) + (n / 2), std::end(SA));
  auto output_SA = std::ranges::subrange(std::begin(SA), std::begin(SA) + (n / 2));
  lms_suffix_direct_sort_dna<char_type, size_type>(S, output_SA, buffer, m, k, num_threads, states);
  SPDLOG_DEBUG("lms_suffix_direct_sort elapsed {}", sw);

  sw = spdlog::stopwatch{};
  auto SA1 = std::ranges::subrange(SA.begin() + 1, SA.begin() + m);
  put_lms_suffix(S, SA, SA1, states);
  SPDLOG_DEBUG("put_lms_suffix elapsed {}", sw);

  sw = spdlog::stopwatch{};
  induced_sort(S, SA, states);
  SPDLOG_DEBUG("induced_sort elapsed {}", sw);
}

template <typename char_type, typename size_type>
void kiss1_suffix_array(
  const vector<char_type> &S,
  vector<size_type> &SA,
  size_type k = 256u,
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
  SA.resize(n / 2);
  SPDLOG_DEBUG("SA.resize(n / 2) elapsed {}", sw);

  sw = spdlog::stopwatch{};
  auto m = get_lms(S, SA, states);
  SPDLOG_DEBUG("get_lms elapsed {}", sw);

  sw = spdlog::stopwatch{};
  lms_suffix_direct_sort<char_type, size_type>(S, SA, m, k);
  SPDLOG_DEBUG("lms_suffix_direct_sort elapsed {}", sw);

  sw = spdlog::stopwatch{};
  SA.resize(n + 1);
  SPDLOG_DEBUG("SA.resize(n + 1) elapsed {}", sw);

  sw = spdlog::stopwatch{};
  auto SA1 = std::ranges::subrange(SA.begin() + 1, SA.begin() + m);
  put_lms_suffix(S, SA, SA1, states);
  SPDLOG_DEBUG("put_lms_suffix elapsed {}", sw);

  sw = spdlog::stopwatch{};
  induced_sort(S, SA, states);
  SPDLOG_DEBUG("induced_sort elapsed {}", sw);
}

} // namespace kiss

} // namespace biovoltron
