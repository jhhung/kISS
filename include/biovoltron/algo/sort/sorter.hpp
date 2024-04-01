#pragma once

#include <biovoltron/utility/istring.hpp>

namespace biovoltron {

template<typename T>
concept SASorter = requires(T t, istring_view ref) {
  t.get_suffix_array_dna(ref);
};

}  // namespace biovoltron
