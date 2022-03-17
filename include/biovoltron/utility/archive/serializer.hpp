#pragma once

#include <cassert>
#include <fstream>
#include <ranges>
#include <vector>

namespace biovoltron {

/**
 * @ingroup utility
 * A C++20 generalized serializer that can serialize the 
 * [continuous_range][continuous_range] which `value_type` is 
 * [trivially_copyable][trivially_copyable], such as scalar type, 
 * simple `struct`, `array` etc. The binary archive file size is 
 * actually the memory occupation in run time plus one `size_t` which 
 * indicates the size of the continuous container.
 * 
 * Furthermore, biovoltron::Serializer can also serialize the 
 * `vector<bool>` and still make the binary archive file size as 
 * small as the original memory occupation. This is hugely different 
 * from [Boost.Serialization][Boost.Serialization] which store the 
 * `vector<bool>` as a `vector<uint8_t>` even in `binary_archive` mode. 
 * 
 * The biovoltron::Serializer can also serialize my other project 
 * biovoltron::DibitVector which is a space-efficient specialization 
 * of the `vector` of two-bits numbers.
 * 
 * Usage
 * ```cpp
 * #include <string>
 * #include <biovoltron/container/xbit_vector.hpp>
 * #include <biovoltron/utility/archive/serializer.hpp>
 * 
 * template <typename R>
 * void test(const R& r, const char* path) {
 *   using biovoltron::Serializer;
 *   {
 *     auto fout = std::ofstream{path, std::ios::binary};
 *     Serializer::save(fout, r);
 *   }
 *   auto r2 = R{};
 *   {
 *     auto fin = std::ifstream{path, std::ios::binary};
 *     Serializer::load(fin, r2);
 *   }
 *   assert(r == r2);
 * }
 * 
 * struct Point {
 *   double x;
 *   double y;
 * };
 * 
 * int main() {
 *   test(std::vector{true, false, false}, "vector_bool.bin");
 *   test(std::vector{7, 5, 16, 8}, "vector_int.bin");
 *   test(std::string{"Exemplar"}, "string.bin");
 *   test(std::vector{Point{0.0, 0.0}, {1.0, 2.0}, {3.0, 4.0}}, "vector_point.bin");
 *   using biovoltron::DibitVector;
 *   test(DibitVector{1, 0, 2, 3, 3, 0, 2}, "DibitVector.bin");
 * }
 * ```
 * 
 * [continuous_range]: https://en.cppreference.com/w/cpp/ranges/contiguous_range
 * [trivially_copyable]: https://en.cppreference.com/w/cpp/types/is_trivially_copyable
 * [Boost.Serialization]: https://www.boost.org/doc/libs/1_76_0/libs/serialization/doc/index.html
 */
struct Serializer {
 private:
  template<std::ranges::random_access_range R>
  constexpr static auto
  get_bytes(const R& r) noexcept {
    if constexpr (requires { r.num_blocks(); })
      return r.num_blocks() * sizeof(typename R::block_type);
    else if constexpr (std::same_as<R, std::vector<bool>>)
      return ((r.size() - 1) / 64 + 1) * sizeof(std::uint64_t);
    else
      return r.size() * sizeof(typename R::value_type);
  }

  template<std::ranges::random_access_range R>
  static auto
  get_data(R& r) noexcept {
    if constexpr (std::same_as<std::remove_const_t<R>, std::vector<bool>>)
      return *reinterpret_cast<std::uint64_t* const*>(&r);
    else
      return r.data();
  }

 public:
  template<std::ranges::random_access_range R, auto BUF_SIZE = 8192>
    requires std::is_trivially_copyable_v<std::ranges::range_value_t<R>>
  static auto
  save(std::ofstream& fout, const R& r) {
    const auto size = r.size();
    if (size == 0)
      return;
    // write size of range
    fout.write(reinterpret_cast<const char*>(&size), sizeof(size));
    const auto bytes = get_bytes(r);
    // write contents
    auto begin = reinterpret_cast<const char*>(get_data(r));
    for (auto batch = 0u; batch < bytes / BUF_SIZE; batch++, begin += BUF_SIZE)
      fout.write(begin, BUF_SIZE);
    // write remains
    if (const auto remains = bytes % BUF_SIZE; remains)
      fout.write(begin, remains);
  }

  template<std::ranges::random_access_range R, auto BUF_SIZE = 8192>
    requires std::is_trivially_copyable_v<std::ranges::range_value_t<R>>
  static auto
  load(std::ifstream& fin, R& r) {
    const auto size = [&] {
      auto size = r.size();
      // read size of results
      fin.read(reinterpret_cast<char*>(&size), sizeof(size));
      return size;
    }();
    if (size == 0)
      return;
    r.resize(size);

    const auto bytes = get_bytes(r);
    // read contents
    auto begin = reinterpret_cast<char*>(get_data(r));
    for (auto batch = 0u; batch < bytes / BUF_SIZE;
         batch++, begin += BUF_SIZE) {
      fin.read(begin, BUF_SIZE);
      assert(static_cast<std::size_t>(fin.gcount()) == BUF_SIZE);
    }
    // read remains
    if (const auto remains = bytes % BUF_SIZE; remains) {
      fin.read(begin, remains);
      assert(static_cast<std::size_t>(fin.gcount()) == remains);
    }
  }
};

}  // namespace biovoltron
