#pragma once

#include <biovoltron/algo/sort/core/psais.hpp>
#include <biovoltron/utility/istring.hpp>
#include <spdlog/spdlog.h>
#include <spdlog/stopwatch.h>

#include <vector>

namespace biovoltron {

template<typename size_type = std::uint32_t>
struct PsaisSorter {
  static auto
  get_sa(istring_view ref, size_t sort_len = istring_view::npos) {
    auto SA = std::vector<size_type>(ref.size() + 1, psais::EMPTY<size_type>);
    auto T = psais::TypeVector(ref.size(), psais::SUFFIX_TYPE::L_TYPE);
    suffix_array(ref, 5u, SA, T);
    return SA;
  }

  static auto
  get_sa(std::string_view ref, size_t sort_len = istring_view::npos) {
    auto SA = std::vector<size_type>(ref.size() + 1, psais::EMPTY<size_type>);
    auto T = psais::TypeVector(ref.size(), psais::SUFFIX_TYPE::L_TYPE);
    suffix_array(ref, 128u, SA, T);
    return SA;
  }

 private:
  static void
  suffix_array(
    const std::ranges::random_access_range auto& S,
    size_type K,
    std::ranges::random_access_range auto& SA,
    auto& T
  ) {
    // 1. get type && bucket
    auto get_type_and_bucket_1_sw = spdlog::stopwatch{};
    psais::get_type<size_type>(S, T);
    auto BA = psais::get_bucket(S, K);

    // 2. put LMS character into SA in any order for each bucket
    psais::put_lms_substr(S, T, BA, SA);
    SPDLOG_DEBUG("GetSuffixType elapsed {}", get_type_and_bucket_1_sw);

    // 3. induce LMS substring
    auto induce_sort_1_sw = spdlog::stopwatch{};
    psais::induce_sort(S, T, BA, SA);
    SPDLOG_DEBUG("InduceSort_1 elapsed {}", induce_sort_1_sw);
    
    BA.clear();
    BA.shrink_to_fit();

    // 4. naming LMS substring
    // |S1| = n1, |SA1| = n1 + 1

    auto name_lms_substr_sw = spdlog::stopwatch{};
    auto [n1, K1] = psais::name_lms_substr<size_type>(S, T, SA);
    SPDLOG_DEBUG("LMSSuffixNaming elapsed {}", name_lms_substr_sw);
    // SPDLOG_DEBUG("InduceSort_1_(include_naming) elapsed {}", induce_sort_1_sw);

    auto S1 = std::ranges::subrange(std::end(SA) - n1, std::end(SA));
    auto SA1 = std::ranges::subrange(std::begin(SA), std::begin(SA) + n1 + 1);
    auto T1 = std::ranges::subrange(std::begin(T), std::begin(T) + n1);

    // 5. recursively solve LMS suffix
    if (K1 < n1) {
      suffix_array(S1, K1, SA1, T1);
    } else {
#pragma omp parallel for
      for (auto i = size_type{}; i < n1; i++) SA1[S1[i] + 1] = i;
      SA1[0] = n1;
    }

    // 6. get type && bucket
    auto get_type_and_bucket_2_sw = spdlog::stopwatch{};
    psais::get_type<size_type>(S, T);
    BA = psais::get_bucket(S, K);
    // SPDLOG_DEBUG("get_type_and_bucket_2 elapsed {}", get_type_and_bucket_2_sw);

    // 7. put LMS character into SA in the order of LMS suffix
    // auto put_lms_suffix_sw = spdlog::stopwatch{};
    psais::put_lms_suffix<size_type>(S, T, BA, S1, SA, SA1);
    // SPDLOG_DEBUG("put_lms_suffix elapsed {}", put_lms_suffix_sw);

    // 8. induce SA from LMS suffix
    // auto induce_sort_2_sw = spdlog::stopwatch{};
    psais::induce_sort(S, T, BA, SA);
    // SPDLOG_DEBUG("induce_sort_2 elapsed {}", induce_sort_2_sw);
    SPDLOG_DEBUG("InduceSort_2 elapsed {}", get_type_and_bucket_2_sw);
  }

};

}  // namespace biovoltron
