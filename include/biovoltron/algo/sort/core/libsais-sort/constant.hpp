#pragma once

#include <limits>
#include <cstdint>
#include <cstddef>

namespace sais {

enum SUFFIX_TYPE {
  L_TYPE   = 0b0,
  S_TYPE   = 0b1,
  LMS_TYPE = 0b10
};
constexpr uint8_t SUFFIX_TYPE_MASK = 0b11;

template <typename T>
constexpr auto EMPTY = std::numeric_limits<T>::max();

constexpr size_t CHAR_SIZE = 256;
constexpr size_t THREAD_CACHE_SIZE = 24 * 1024;

}
