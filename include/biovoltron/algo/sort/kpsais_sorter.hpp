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
    SA.resize(n1 + 1, psais::EMPTY<size_type>);

    // 3. place lms index
    auto buf = std::ranges::subrange(std::begin(SA) + 1, std::end(SA));
    psais::put_lms_suffix_left_shift<size_type>(T, buf);
    SA[0] = n;
    SPDLOG_DEBUG("GetLMSPositions elapsed {}", sw);

    // 4. sort lms suffix in sort_len order
    auto reverse_bits = [](uint16_t x) {
      x = ((x >> 1) & 0x5555) | ((x & 0x5555) << 1);
      x = ((x >> 2) & 0x3333) | ((x & 0x3333) << 2);
      x = ((x >> 4) & 0x0f0f) | ((x & 0x0f0f) << 4);
      x = ((x >> 8) & 0x00ff) | ((x & 0x00ff) << 8);
      return x;
    };
    std::vector<uint16_t> reverse_bits_table;
    for (int i = 0; i < (1 << 16); i++)
      reverse_bits_table.emplace_back(reverse_bits(i));
    sw = spdlog::stopwatch{};
    std::stable_sort(
      std::execution::par_unseq, std::begin(SA), std::end(SA),
      [&ref, sort_len, n, &reverse_bits_table](size_type i, size_type j) {
        //  return ref.substr(i, sort_len) < ref.substr(j, sort_len);
        auto tmp = size_type{};
        for (; tmp + 16 <= sort_len && i + 16 <= n && j + 16 <= n;
             tmp += 16, i += 16, j += 16) {
          __m128i a
            = _mm_loadu_si128(reinterpret_cast<const __m128i*>(&ref[i]));
          __m128i b
            = _mm_loadu_si128(reinterpret_cast<const __m128i*>(&ref[j]));
          const uint16_t less = _mm_movemask_epi8(_mm_cmplt_epi8(a, b));
          const uint16_t greater = _mm_movemask_epi8(_mm_cmpgt_epi8(a, b));
          if (less || greater) {
            return reverse_bits_table[less] > reverse_bits_table[greater];
          }
        }
        for (; tmp + 1 <= sort_len && i < n && j < n; tmp++, i++, j++) {
          if (ref[i] < ref[j])
            return true;
          if (ref[i] > ref[j])
            return false;
        }
        if (tmp >= sort_len) {
          return i < j;
        }
        if (i == n)
          return true;
        else
          return false;
      });
    SPDLOG_DEBUG("parallelkOrderedLMSSuffixSort elapsed {}", sw);
    // 5. get bucket
    sw = spdlog::stopwatch{};
    auto BA = psais::get_bucket(ref, K);
    SPDLOG_DEBUG("getCharacterBucket elapsed {}", sw);

    // 6. induce SA
    sw = spdlog::stopwatch{};
    SA.resize(n + 1, psais::EMPTY<size_type>);
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
