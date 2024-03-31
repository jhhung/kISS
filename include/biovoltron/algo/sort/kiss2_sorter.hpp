#pragma once
// TODO: spdlog
// TODO: unify include format
#include <biovoltron/algo/sort/kiss2_core.hpp>
#include <vector>

namespace biovoltron {

template<typename size_type = std::uint32_t>
struct KISS2Sorter {
  static auto
  get_suffix_array_dna(const std::ranges::random_access_range auto& ref,
        size_type k = 256u,
        const size_t num_threads = std::thread::hardware_concurrency()) {
    auto S = kiss::vector<uint8_t>{ref.begin(), ref.end()};
    auto SA = kiss::vector<size_type>{};
    kiss::kiss2_suffix_array_dna(S, SA, k, num_threads);
    return std::vector<size_type>(std::make_move_iterator(SA.begin()), std::make_move_iterator(SA.end()));
  }

  static auto
  get_suffix_array(const std::ranges::random_access_range auto& ref,
        size_type k = 256u,
        const size_t num_threads = std::thread::hardware_concurrency()) {
    auto S = kiss::vector<uint8_t>{ref.begin(), ref.end()};
    auto SA = kiss::vector<size_type>{};
    kiss::kiss2_suffix_array(S, SA, k, num_threads);
    return std::vector<size_type>(std::make_move_iterator(SA.begin()), std::make_move_iterator(SA.end()));
  }
};

}
