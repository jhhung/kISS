#pragma once

#include <biovoltron/utility/istring.hpp>
#include <execution>
#include <numeric>

namespace biovoltron {

struct Sorter {
  template<typename size_type = std::uint32_t>
  static auto
  get_sa(istring_view ref, std::size_t sort_len = istring_view::npos) {
    auto sa = std::vector<size_type>(ref.size() + 1);
    std::iota(sa.begin(), sa.end(), 0);
    std::stable_sort(std::execution::par_unseq, sa.begin(), sa.end(),
                     [ref, sort_len](auto i, auto j) {
                       return ref.substr(i, sort_len) < ref.substr(j, sort_len);
                     });
    return sa;
  }
};

}  // namespace biovoltron