#pragma once
#include <biovoltron/algo/sort/kiss1_core.hpp>
#include <vector>

namespace biovoltron {

template<typename size_type = std::uint32_t>
struct KISS1Sorter {
  using SA_t = kiss::vector<size_type>;
  static auto
  get_suffix_array_dna(const std::ranges::random_access_range auto& ref,
        size_type k = 256u,
        const size_t num_threads = std::thread::hardware_concurrency()) {
    auto S = kiss::vector<uint8_t>{ref.begin(), ref.end()};
    auto SA = kiss::vector<size_type>{};
    kiss::kiss1_suffix_array_dna(S, SA, k, num_threads);
    return SA;
  }

  static auto
  get_suffix_array_dna(const kiss::vector<uint8_t>& S, size_type k = 256u,
        const size_t num_threads = std::thread::hardware_concurrency()) {
    auto SA = kiss::vector<size_type>{};
    kiss::kiss1_suffix_array_dna<uint8_t, size_type>(S, SA, k, num_threads);
    return SA;
  }

  static auto
  get_suffix_array(const std::ranges::random_access_range auto& ref,
        size_type k = 256u,
        const size_t num_threads = std::thread::hardware_concurrency()) {
    auto S = kiss::vector<uint8_t>{ref.begin(), ref.end()};
    auto SA = kiss::vector<size_type>{};
    kiss::kiss1_suffix_array(S, SA, k, num_threads);
    return SA;
  }

  static auto
  get_suffix_array(const kiss::vector<uint8_t>& S, size_type k = 256u,
        const size_t num_threads = std::thread::hardware_concurrency()) {
    auto SA = kiss::vector<size_type>{};
    kiss::kiss1_suffix_array<uint8_t, size_type>(S, SA, k, num_threads);
    return SA;
  }

  static auto
  prepare_aligned_ref(const std::ranges::random_access_range auto& ref) {
    return kiss::vector<uint8_t>{ref.begin(), ref.end()};
  }
};

}
