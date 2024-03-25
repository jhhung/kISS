#pragma once

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
    // suffix_array(ref, 5u, SA, sort_len);
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
               std::ranges::random_access_range auto& SA, size_type sort_len) {
    size_type n = ref.size();

    // 1. Type Array Generation: processes the reference string S and generates
    // a Boolean array T (referred to as the "type array") with a length of |S|.
    // This array characterizes the type of each suffix in the input string.
    auto sw = spdlog::stopwatch{};
    auto T = psais::TypeVector(n, psais::SUFFIX_TYPE::L_TYPE);
    psais::get_type<size_type>(ref, T);
    SPDLOG_DEBUG("classifyLSType elapsed {}", sw);

    // 2. prepare lms array
    sw = spdlog::stopwatch{};
    auto n1 = psais::num_lms<size_type>(T);

    auto sw3 = spdlog::stopwatch{};
    SA.reserve(n + 1);
    SA.resize(n1 + 1, psais::EMPTY<size_type>);
    SPDLOG_DEBUG("First reserve & resize {}", sw3);

    // 3. place lms index
    auto lms_without_sentinel
      = std::ranges::subrange(std::begin(SA), std::end(SA) - 1);
    psais::put_lms_suffix_left_shift<size_type>(T, lms_without_sentinel);
    SA.back() = n;
    SPDLOG_DEBUG("GetLMSPositions elapsed {}", sw);

    // 4. compress the string + get LMS positions
    // FIXME: may suffer from errors if n2 * 4 + 2 does not fit into
    // size_type
    sw = spdlog::stopwatch{};
    // Consecutive l-mers will be encoded into one integer.
    auto sw2 = spdlog::stopwatch{};
    size_type l = 8 * sizeof(size_type) / std::bit_width(K - 1);
    // The length of the encoded string is precalculated to estimate memory
    // usage precisely.
    auto n2 = kiss::get_encoded_reference_length<size_type>(SA, l);
    SPDLOG_DEBUG("Calculate encode length {}", sw2);
    sw2 = spdlog::stopwatch{};
    SA.resize(std::max(n2 * 3 + 1, n + 1));
    SPDLOG_DEBUG("Resize {}", sw2);
    sw2 = spdlog::stopwatch{};
    // Identify the range for lms (LMS), rank (Inverse suffix array), buf
    // (buffer), stating_position and valid_position
    auto lms = std::ranges::subrange(std::begin(SA), std::begin(SA) + (n1 + 1));
    auto rank = std::ranges::subrange(std::begin(SA) + (n2 + 1),
                                      std::begin(SA) + (n2 * 2 + 1));
    auto buf = std::ranges::subrange(std::begin(SA) + (n2 * 2 + 1),
                                     std::begin(SA) + (n2 * 2 + 1) + n1);
    // auto starting_position = std::ranges::subrange(
    //   std::begin(SA) + (n2 * 3 + 1), std::begin(SA) + (n2 * 4 + 2));
    auto valid_position = psais::TypeVector(n2, 0);
    kiss::encode_reference<size_type>(ref, lms, rank, buf,
                                     valid_position, K, l);
    SPDLOG_DEBUG("Encode reference elapsed {}", sw2);

    // 5. prefix doubling & place back
    // prefix doubling
    auto sa = std::ranges::subrange(std::begin(SA) + (n2 * 2 + 1),
                                    std::begin(SA) + (n2 * 3 + 1));
    buf = std::ranges::subrange(std::begin(SA) + 1, std::begin(SA) + (n2 + 1));
    kiss::prefix_doubling(sa, rank, buf, sort_len);
    // place back
    sw2 = spdlog::stopwatch{};
    auto sorted_lms = buf;
    psais::put_lms_suffix_left_shift<size_type>(T, rank);
    // rank[n1] = n;
    lms = rank;
    kiss::place_back_lms<size_type>(sa, sorted_lms, lms,
                                    valid_position, l, n1 + 1, n);
    SPDLOG_DEBUG("Place back LMS {}", sw2);
    // place sorted_lms to the beginning of SA
    if (std::begin(sorted_lms) - (std::begin(SA) + 1))
      std::copy(std::begin(sorted_lms),
                std::begin(sorted_lms) + n1, std::begin(SA) + 1);
    SA[0] = n;
    SPDLOG_DEBUG("parallelkOrderedLMSSuffixSort elapsed {}", sw);

    // 6. get bucket
    sw = spdlog::stopwatch{};
    auto BA = psais::get_bucket(ref, K);
    SPDLOG_DEBUG("getCharacterBucket elapsed {}", sw);

    // 7. induce SA
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
