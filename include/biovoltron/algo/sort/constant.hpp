#pragma once

#include <algorithm>
#include <limits>
#include <cstdint>
#include <cstddef>

namespace biovoltron {

namespace kiss {

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

constexpr size_t TYPEVECTOR_ELEMENT_SIZE = 8;

constexpr size_t KISS1_SPLIT_SORT_PREFIX_SIZE = 10;
constexpr size_t KISS1_SPLIT_SORT_BUCKET_SIZE = (1 << (2 * KISS1_SPLIT_SORT_PREFIX_SIZE));
constexpr size_t KISS1_SPLIT_SORT_STRIDE_DNA = 125;
constexpr size_t KISS1_SPLIT_SORT_STRIDE = 32;

constexpr size_t KISS2_RADIX_SORT_BUCKET_SIZE = (1 << 16);
constexpr size_t KISS2_RADIX_SORT_MASK = ((1 << 16) - 1);
constexpr int KISS2_HEAD_MASK = 2;
constexpr int KISS2_NEW_HEAD_MASK = 1;

constexpr size_t THREADSTATE_MAX_BUFFER_SIZE = std::max(KISS1_SPLIT_SORT_BUCKET_SIZE, KISS2_RADIX_SORT_BUCKET_SIZE);

} // namespace kiss

} // namespace biovoltron