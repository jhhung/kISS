#pragma once

#include <biovoltron/algo/sort/core/psais.hpp>
#include <biovoltron/utility/istring.hpp>

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
    psais::get_type<size_type>(S, T);
    auto BA = psais::get_bucket(S, K);

    // 2. put LMS character into SA in any order for each bucket
    psais::put_lms_substr(S, T, BA, SA);

    // 3. induce LMS substring
    psais::induce_sort(S, T, BA, SA);
    BA.clear();
    BA.shrink_to_fit();

    // 4. naming LMS substring
    // |S1| = n1, |SA1| = n1 + 1
    auto [n1, K1] = psais::name_lms_substr<size_type>(S, T, SA);
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
    psais::get_type<size_type>(S, T);
    BA = psais::get_bucket(S, K);

    // 7. put LMS character into SA in the order of LMS suffix
    psais::put_lms_suffix<size_type>(S, T, BA, S1, SA, SA1);

    // 8. induce SA from LMS suffix
    psais::induce_sort(S, T, BA, SA);
  }

};

}  // namespace biovoltron
