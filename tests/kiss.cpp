#include <biovoltron/algo/sort/kiss1_sorter.hpp>
#include <biovoltron/algo/sort/kiss_new_sorter.hpp>
#include <biovoltron/utility/istring.hpp>
#include <experimental/random>
#include <catch_amalgamated.hpp>
#include <chrono>
#include <algorithm>

using namespace biovoltron;

TEST_CASE("kISS-1 Sorter") {
  auto gen_dna_seq = [](int len) {
    auto seq = std::string{};
    while (len--)
      seq += "atgc"[std::experimental::randint(0, 3)];
    return seq;
  };

  int len = std::experimental::randint(100'000, 200'000);
  auto k = size_t{256};
  const auto seq = gen_dna_seq(len);
  const auto seq_sv = std::string_view{seq};
  const auto ref = Codec::to_istring(seq);

  auto sa = KISS1Sorter<>::get_sa(ref, k);
#pragma omp parallel for
  for (int i = 1; i < sa.size(); i++)
    REQUIRE(seq_sv.substr(sa[i - 1], k) <= seq_sv.substr(sa[i], k));
}

TEST_CASE("kISS-1 large testcase") {
  auto gen_dna_seq = [](int len) {
    auto seq = std::string{};
    while (len--)
      seq += "atgc"[std::experimental::randint(0, 3)];
    return seq;
  };

  int len = std::experimental::randint(10'000'000, 20'000'000);
  auto k = size_t{256};
  const auto seq = gen_dna_seq(len);
  const auto seq_sv = std::string_view{seq};
  const auto ref = Codec::to_istring(seq);

  auto sa = KISS1Sorter<>::get_sa(ref, k);
#pragma omp parallel for
  for (int i = 1; i < sa.size(); i++)
    REQUIRE(seq_sv.substr(sa[i - 1], k) <= seq_sv.substr(sa[i], k));
}

TEST_CASE("kISS-2 Sorter") {
  auto gen_dna_seq = [](int len) {
    auto seq = std::string{};
    while (len--)
      seq += "atgc"[std::experimental::randint(0, 3)];
    return seq;
  };

  int len = std::experimental::randint(100'000, 200'000);
  auto k = size_t{256};
  const auto seq = gen_dna_seq(len);
  const auto seq_sv = std::string_view{seq};
  const auto ref = Codec::to_istring(seq);

  auto sa = KissNewSorter<>::get_sa(ref, k);
#pragma omp parallel for
  for (int i = 1; i < sa.size(); i++)
    REQUIRE(seq_sv.substr(sa[i - 1], k) <= seq_sv.substr(sa[i], k));
}

TEST_CASE("kISS-2 Sorter large testcase") {
  auto gen_dna_seq = [](int len) {
    auto seq = std::string{};
    while (len--)
      seq += "atgc"[std::experimental::randint(0, 3)];
    return seq;
  };

  int len = std::experimental::randint(10'000'000, 20'000'000);
  auto k = size_t{256};
  const auto seq = gen_dna_seq(len);
  const auto seq_sv = std::string_view{seq};
  const auto ref = Codec::to_istring(seq);

  auto sa = KissNewSorter<>::get_sa(ref, k);
#pragma omp parallel for
  for (int i = 1; i < sa.size(); i++)
    REQUIRE(seq_sv.substr(sa[i - 1], k) <= seq_sv.substr(sa[i], k));
}
