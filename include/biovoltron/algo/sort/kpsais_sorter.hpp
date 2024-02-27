#pragma once

#include <immintrin.h>
#include <spdlog/spdlog.h>
#include <spdlog/stopwatch.h>

#include <biovoltron/algo/sort/core/psais.hpp>
#include <biovoltron/algo/sort/core/sorter.hpp>
#include <biovoltron/container/xbit_vector.hpp>
#include <biovoltron/utility/istring.hpp>
#include <execution>
#include <thread>

namespace biovoltron {

template<typename size_type = std::uint32_t>
struct KPsaisSorter {
  static auto
  get_sa(istring_view ref, size_t sort_len = 256u) {
    auto SA = std::vector<size_type>{};
    suffix_array(ref, 4u, SA, sort_len);
    return SA;
  }

  static auto
  get_sa(std::string_view ref, size_t sort_len = 256u) {
    auto SA = std::vector<size_type>{};
    suffix_array(ref, 128u, SA, sort_len);
    return SA;
  }

 private:
  static void
  suffix_array(const std::ranges::random_access_range auto& ref, size_type K,
               std::ranges::random_access_range auto& SA, size_t sort_len) {
    auto n = ref.size();

    // 1. get type
    auto sw = spdlog::stopwatch{};
    auto T = psais::TypeVector(n, psais::SUFFIX_TYPE::L_TYPE);
    psais::get_type<size_type>(ref, T);
    SPDLOG_DEBUG("classifyLSType elapsed {}", sw);

    // 2. prepare lms array
    sw = spdlog::stopwatch{};
    auto n1 = psais::num_lms<size_type>(T);
    SA.reserve(n + 1);
    SA.resize(2 * (n1 + 1), psais::EMPTY<size_type>);

    // 3. place lms index
    auto buf_SA = std::ranges::subrange(std::begin(SA) + (n1 + 2), std::end(SA));
    psais::put_lms_suffix_left_shift<size_type>(T, buf_SA);
    SA[n1 + 1] = n;
    buf_SA = std::ranges::subrange(std::begin(SA) + (n1 + 1), std::end(SA));
    SPDLOG_DEBUG("GetLMSPositions elapsed {}", sw);

    // 4. sort lms suffix in sort_len order
    sw = spdlog::stopwatch{};
    auto output_SA = std::ranges::subrange(std::begin(SA), std::begin(SA) + (n1 + 1));
    const auto prefix_size = 10;
    psais::parallel_k_ordered_sort<size_type>(ref, buf_SA, output_SA, sort_len, prefix_size);
    // auto reverse_packed = psais::get_reverse_packed<size_type>(ref);
    // std::stable_sort(std::execution::par_unseq,
    //   std::begin(SA), std::end(SA),
    //   [&ref, sort_len, n, &reverse_packed](size_type i, size_type j) {
    //     const auto stride = 125;
    //     auto sorted_len = size_type{};
    //     for (; sorted_len + stride <= sort_len && i + stride <= n && j + stride <= n;
    //          sorted_len += stride, i += stride, j += stride) {
    //       auto a = psais::load_125<size_type>(reverse_packed, n, i);
    //       auto b = psais::load_125<size_type>(reverse_packed, n, j);
    //       auto eq = _mm256_cmpeq_epi8(a, b);
    //       auto neq_mask = ~((uint32_t)_mm256_movemask_epi8(eq));

    //       if (neq_mask != 0) {
    //         auto msb_mask = (1UL << (31 - _lzcnt_u32(neq_mask)));
    //         auto gt = _mm256_max_epu8(a, b);
    //         gt = _mm256_cmpeq_epi8(gt, a);
    //         auto gt_mask = (uint32_t)_mm256_movemask_epi8(gt);
    //         if ((msb_mask & gt_mask) > 0) {
    //             return false;
    //         } else {
    //           return true;
    //         }
    //       }
    //     }
    //     for (; sorted_len + 1 <= sort_len && i < n && j < n; sorted_len++, i++, j++) {
    //       if (ref[i] < ref[j])
    //         return true;
    //       if (ref[i] > ref[j])
    //         return false;
    //     }
    //     if (sorted_len >= sort_len) {
    //       return i < j;
    //     }
    //     if (i == n)
    //       return true;
    //     else
    //       return false;
      // });
    SPDLOG_DEBUG("parallelkOrderedLMSSuffixSort elapsed {}", sw);
    // 5. get bucket
    sw = spdlog::stopwatch{};
    auto BA = psais::get_bucket(ref, K);
    SPDLOG_DEBUG("getCharacterBucket elapsed {}", sw);

    // 6. induce SA
    sw = spdlog::stopwatch{};
    SA.resize(n + 1, psais::EMPTY<size_type>);
    // set the rest of SA to default value
    auto remain_len = (n - n1);
#pragma omp parallel for num_threads(std::thread::hardware_concurrency())
    for (auto i = size_type{}; i < remain_len; i++) {
      SA[i + (n1 + 1)] = psais::EMPTY<size_type>;
    }
    auto ref1 = std::ranges::subrange(std::end(SA) - n1, std::end(SA));
    auto SA1 = std::ranges::subrange(std::begin(SA), std::begin(SA) + n1 + 1);
    psais::put_lms_suffix_right_shift<size_type>(ref, T, BA, ref1, SA, SA1);
    SPDLOG_DEBUG("putLMSSuffixes elapsed {}", sw);
    sw = spdlog::stopwatch{};
    psais::induce_sort(ref, T, BA, SA);
    SPDLOG_DEBUG("InduceSort elapsed {}", sw);
  }
};

}  // namespace biovoltron
