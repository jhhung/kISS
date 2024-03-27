#pragma once

#include "constant.hpp"

#include <vector>

template<typename T, std::size_t Alignment>
class AlignedAllocator {
 public:
  using value_type = T;
  using pointer = T*;
  using const_pointer = const T*;
  using reference = T&;
  using const_reference = const T&;
  using size_type = std::size_t;
  using difference_type = std::ptrdiff_t;

  template<typename U>
  struct rebind {
    using other = AlignedAllocator<U, Alignment>;
  };

  AlignedAllocator() noexcept {}

  template<typename U>
  AlignedAllocator(const AlignedAllocator<U, Alignment>&) noexcept {}

  pointer allocate(size_type n) {
    if (auto p = aligned_alloc(Alignment, n * sizeof(T)))
      return static_cast<pointer>(p);
    throw std::bad_alloc();
  }

  void deallocate(pointer p, size_type) noexcept {
    std::free(p);
  }

  friend bool operator==(const AlignedAllocator&, const AlignedAllocator&) noexcept {
    return true;
  }

  friend bool operator!=(const AlignedAllocator& a, const AlignedAllocator& b) noexcept {
    return !(a == b);
  }
};

namespace sais {

template <typename T>
using vector = std::vector<T, AlignedAllocator<T, 64>>;

template <typename size_type>
struct ThreadCache {
  size_type symbol, index;
};

template <typename size_type>
struct alignas(64) ThreadState {
  size_type beg;
  size_type len;
  size_type m;
  sais::vector<size_type> bucket;
  sais::vector<ThreadCache<size_type>> cache;

  ThreadState() : bucket(5 * CHAR_SIZE), cache(THREAD_CACHE_SIZE) {}
};

}
