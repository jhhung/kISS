#pragma once

#include <bits/stdc++.h>
#include <spdlog/spdlog.h>
#include <spdlog/stopwatch.h>

#include <biovoltron/algo/sort/core/kiss_new.hpp>
#include <biovoltron/algo/sort/core/psais.hpp>
#include <biovoltron/algo/sort/core/sorter.hpp>
#include <biovoltron/container/xbit_vector.hpp>
#include <biovoltron/utility/istring.hpp>
#include <bit>
#include <execution>
#include <thread>

namespace biovoltron {

template<typename size_type = std::uint32_t>
struct KissNewSorter {
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
               std::ranges::random_access_range auto& SA, size_type sort_len) {
    auto n = (size_type)ref.size();

    // 1. get type
    auto GetLMSPosition_sw = spdlog::stopwatch{};
    auto T = psais::TypeVector(n, psais::SUFFIX_TYPE::L_TYPE);
    psais::get_type<size_type>(ref, T);

    // 2. prepare lms array
    auto n1 = psais::num_lms<size_type>(T);
    auto compress_block_length
      = (size_type)(8 * sizeof(size_type) / std::bit_width(K - 1));
    auto n2 = kiss::len_compressed_string<size_type>(T, compress_block_length);
    SA.reserve(std::max(n2 * 4 + 2, n + 1));
    SA.resize(n1 + 1, psais::EMPTY<size_type>);

    // 3. place lms index
    auto buf = std::ranges::subrange(std::begin(SA), std::end(SA) - 1);
    psais::put_lms_suffix_left_shift<size_type>(T, buf);
    SA.back() = n;
    SPDLOG_DEBUG("GetSuffixType elapsed {}", GetLMSPosition_sw);

    // 4. compress the string + get LMS positions
    // FIXME: may suffer from errors if n2 * 4 + 2 does not fit into
    // size_type
    auto new_sw = spdlog::stopwatch{};
    auto cp1 = spdlog::stopwatch{};
    SA.resize(std::max(n2 * 4 + 2, n + 1));
    auto lms = std::ranges::subrange(std::begin(SA), std::begin(SA) + (n1 + 1));
    auto rank = std::ranges::subrange(std::begin(SA) + (n2 + 1),
                                      std::begin(SA) + (n2 * 2 + 1));
    buf = std::ranges::subrange(std::begin(SA) + (n2 * 2 + 1),
                                std::begin(SA) + (n2 * 2 + 1) + n1);
    auto starting_position = std::ranges::subrange(
      std::begin(SA) + (n2 * 3 + 1), std::begin(SA) + (n2 * 4 + 2));
    auto valid_position = psais::TypeVector(n2, 0);
    kiss::compress_string<size_type>(ref, lms, rank, buf, starting_position,
                                     valid_position, K, compress_block_length);
    SPDLOG_DEBUG("Compress string elapsed {}", cp1);

    // 5. prefix doubling & place back
    auto sa = std::ranges::subrange(std::begin(SA) + (n2 * 2 + 1),
                                    std::begin(SA) + (n2 * 3 + 1));
    buf = std::ranges::subrange(std::begin(SA) + 1, std::begin(SA) + (n2 + 1));
    kiss::prefix_doubling(sa, rank, buf, starting_position, sort_len);
    // place back
    auto cp3 = spdlog::stopwatch{};
    auto sorted_lms = buf;
    buf = rank;
    kiss::place_back_lms<size_type>(T, sa, sorted_lms, buf, starting_position,
                                    valid_position);
    if (std::begin(sorted_lms) - (std::begin(SA) + 1))
      std::copy(std::execution::par, std::begin(sorted_lms),
                std::begin(sorted_lms) + n1, std::begin(SA) + 1);
    SA[0] = n;
    SPDLOG_DEBUG("Place back elapsed {}", cp3);
    SPDLOG_DEBUG("Prefix Doubling elapsed {}", new_sw);

    // 6. get bucket
    auto InduceSort_sw = spdlog::stopwatch{};
    auto BA = psais::get_bucket(ref, K);

    // 7. induce SA
    SA.resize(n1 + 1);
    SA.shrink_to_fit();
    SA.resize(n + 1, psais::EMPTY<size_type>);
    auto ref1 = std::ranges::subrange(std::end(SA) - n1, std::end(SA));
    auto SA1 = std::ranges::subrange(std::begin(SA), std::begin(SA) + n1 + 1);
    psais::put_lms_suffix_right_shift<size_type>(ref, T, BA, ref1, SA, SA1);
    psais::induce_sort(ref, T, BA, SA);
    SPDLOG_DEBUG("InduceSort elapsed {}", InduceSort_sw);
  }
};

}  // namespace biovoltron
