#pragma once

#include <bits/stdc++.h>
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
    suffix_array(ref, 5u, SA, sort_len);
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
    auto GetLMSPosition_sw = spdlog::stopwatch{};
    auto T = psais::TypeVector(n, psais::SUFFIX_TYPE::L_TYPE);
    psais::get_type<size_type>(ref, T);

    // 2. prepare lms array
    auto n1 = psais::num_lms<size_type>(T);
    SA.reserve(n + 1);
    SA.resize(n1 + 1, psais::EMPTY<size_type>);

    // 3. place lms index
    auto buf = std::ranges::subrange(std::begin(SA) + 1, std::end(SA));
    psais::put_lms_suffix_left_shift<size_type>(T, buf);
    SA[0] = n;
    SPDLOG_DEBUG("GetSuffixType elapsed {}", GetLMSPosition_sw);

    // 4. sort lms suffix in sort_len order
    auto StableSort_sw = spdlog::stopwatch{};
    std::stable_sort(std::execution::par_unseq, std::begin(SA), std::end(SA),
                     [ref, sort_len](size_type i, size_type j) {
                       return ref.substr(i, sort_len) < ref.substr(j, sort_len);
                     });
    SPDLOG_DEBUG("StableSort elapsed {}", StableSort_sw);
    for (auto i : SA) std::cout << i << ' ';
    std::cout << std::endl;

    // 5. get bucket
    auto InduceSort_sw = spdlog::stopwatch{};
    auto BA = psais::get_bucket(ref, K);

    // 6. induce SA
    SA.resize(n + 1, psais::EMPTY<size_type>);
    auto ref1 = std::ranges::subrange(std::end(SA) - n1, std::end(SA));
    auto SA1 = std::ranges::subrange(std::begin(SA), std::begin(SA) + n1 + 1);
    psais::put_lms_suffix_right_shift<size_type>(ref, T, BA, ref1, SA, SA1);
    psais::induce_sort(ref, T, BA, SA);
    SPDLOG_DEBUG("InduceSort elapsed {}", InduceSort_sw);
  }
};

}  // namespace biovoltron
