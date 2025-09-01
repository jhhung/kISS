#pragma once

#include <biovoltron/algo/sort/utils.hpp>
#include <biovoltron/algo/sort/constant.hpp>
#include <biovoltron/container/xbit_vector.hpp>

#include <immintrin.h>
#include <omp.h>
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
    // The required size needs to be a multiple of Alignment; otherwise, round up.
    auto required_size = (n * sizeof(T) + (Alignment - 1)) & ~(Alignment - 1);
    if (auto p = aligned_alloc(Alignment, required_size))
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

namespace biovoltron {

namespace kiss {

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
  vector<size_type> bucket;
  vector<ThreadCache<size_type>> cache;
  vector<size_type> buffer;

  ThreadState() : bucket(5 * CHAR_SIZE), cache(THREAD_CACHE_SIZE), buffer(THREADSTATE_MAX_BUFFER_SIZE) {}
};

/*
This implementation refers to Daniel Liu's "RevPacked" struct
of the simple-saca project:
https://github.com/Daniel-Liu-c0deb0t/simple-saca/
*/
// TODO: check if aligned(64) can be faster
template <typename char_type, typename size_type>
struct alignas(64) PackedDNAString {
  size_type n_padded;
  biovoltron::DibitVector<> packed_dna;

  PackedDNAString(
    const std::ranges::random_access_range auto &S,
    vector<ThreadState<size_type>> &states,
    size_t num_threads
  ) {
    auto n = S.size();
    n_padded = (n + 16 + 3) & ~(size_type(0b11));

    packed_dna.resize(n_padded, 0);

#pragma omp parallel num_threads(num_threads)
    {
      auto tid = omp_get_thread_num();
      auto &state = states[tid];
      distribute_workload_packed_DNA<size_type>(n, state);

      auto block_ptr = packed_dna.data() + ((n_padded - state.beg - 1) / 4);
      auto i = state.beg;

      for (; i + 3 < state.beg + state.len; i += 4, block_ptr--) {
        *block_ptr = (S[i] << 6) | (S[i + 1] << 4) | (S[i + 2] << 2) | (S[i + 3]);
      }

      while (i < state.beg + state.len) {
        packed_dna[n_padded - i - 1] = S[i];
        i++;
      }
    }
  }

  /*
  This function refers to the "load_125" function
  from the reference implementation.
  */
  auto load_prefix_length_125 (size_type idx) {
    idx = n_padded - idx - 128;
    auto i = (idx + 3) / 4;
    auto j = (idx + 3) % 4;
    auto val = _mm256_loadu_si256(reinterpret_cast<const __m256i*>(packed_dna.data() + i));

    auto left_shift = _mm256_set1_epi64x(((3 - j) * 2));
    auto hi = _mm256_sllv_epi64(val, left_shift);
    auto right_shift = _mm256_set1_epi64x(((32 - (3 - j)) * 2));
    auto lo = _mm256_srlv_epi64(_mm256_permute4x64_epi64(val, 0b10010011), right_shift);

    auto mask = _mm256_set_epi8(
      -1,
      -1,
      -1,
      -1,
      -1,
      -1,
      -1,
      -1,
      -1,
      -1,
      -1,
      -1,
      -1,
      -1,
      -1,
      -1,
      -1,
      -1,
      -1,
      -1,
      -1,
      -1,
      -1,
      -1,
      -1,
      -1,
      -1,
      -1,
      -1,
      -1,
      -1,
      0b11000000
    );

   return _mm256_and_si256(_mm256_or_si256(hi, lo), mask);
  }

  /*
  This function refers to the "load_k" function
  from the reference implementation.
  */
  auto load_prefix_length_less_than_16(
    size_type idx,
    size_type prefix_size
  ) {
   idx = n_padded - idx - 16;
   auto i = (idx + 3) / 4;
   auto j = (idx + 3) % 4;
   auto value = *(reinterpret_cast<const uint32_t*>(packed_dna.data() + i));
   return (value << ((3 - j) * 2)) >> ((16 - prefix_size) * 2);
  }
};

using TypeVector
    = biovoltron::detail::XbitVector<1, std::uint8_t, std::allocator<std::uint8_t>>;

} // namespace kiss

} // namespace biovoltron
