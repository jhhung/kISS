#pragma once

#include <spdlog/spdlog.h>

#include <biovoltron/algo/sort/sorter.hpp>
#include <biovoltron/algo/sort/kiss1_sorter.hpp>
#include <biovoltron/container/xbit_vector.hpp>
#include <biovoltron/utility/archive/serializer.hpp>
#include <biovoltron/utility/istring.hpp>
#include <chrono>
#include <execution>
#include <queue>
#include <span>
#include <thread>
#include <utility>

namespace biovoltron {
using namespace std::chrono;

/**
 * @ingroup align
 * @brief
 * A specialization of FM-Index for genome sequence
 * which implemented using C++20.
 *
 * The memory usage for `FMIndex` in run time is affect by above
 * parameters. Take a `3.1Gb` human genome for example, the
 * memory occupation can be calculate as follwing:
 * - bwt string: `3.1Gb / 4 = 0.775 Gb`.
 * - hierarchical occurrence table: L1 occ occupy fixed `3.1Gb * 16 / 256 =
 * 0.194Gb` plus L2 occ occupy `3.1Gb * 4 / occ_intv(16) = 0.775 Gb`.
 * - suffix array: `3.1Gb * 4 / sa_intv(1) = 12Gb`. Noticed that
 *   in default mode we dont sampling suffix value since this can
 *   reduce frequently memory allocation and intense computation
 *   when occurs massive query.
 * - lookup table: `4^lookup_len(14) * 4 / 1024^3 = 1Gb`.
 *
 * So the default total memory occupation for `3.1Gb` human genome
 * is `0.775 Gb + 0.194Gb + 0.775Gb + 12Gb + 1Gb = 14.744Gb`.
 *
 * Example
 * ```cpp
 * #include <iostream>
 * #include <ranges>
 * #include <biovoltron/algo/align/exact_match/fm_index.hpp>
 * #include <biovoltron/file_io/fasta.hpp>
 *
 * using namespace biovoltron;
 *
 * int main(int argc, char** argv) {
 *   if (argc != 2) {
 *     std::cout << "argv[1] should be your fasta file.\n";
 *     return 0;
 *   }
 *
 *   auto fname = std::string(argv[1]);
 *   {
 *     // input
 *     auto ref = istring{};
 *     auto fin = std::ifstream(fname);
 *     for (const auto& r : std::ranges::istream_view<FastaRecord<true>>(fin))
 *       ref += r.seq;
 *
 *     // fm-index only support 'ACTG', so we need to change 'N' to 'ACGT'
 *     std::ranges::transform(ref, ref.begin(), [](auto& c) { return c % 4; });
 *
 *     // build
 *     auto fmi = FMIndex{};
 *     fmi.build(ref);
 *
 *     // save
 *     auto fout = std::ofstream(fname + ".fmi", std::ios::binary);
 *     fmi.save(fout);
 *     std::cout << fname + ".fmi saved.\n";
 *   }
 *
 *   // load
 *   auto fmi = FMIndex{};
 *   {
 *     auto fin = std::ifstream(fname + ".fmi", std::ios::binary);
 *     fmi.load(fin);
 *     std::cout << fname + ".fmi loaded.\n";
 *   }
 *
 *   // query
 *   for (auto seed = istring{}; std::cin >> seed;) {
 *     const auto [beg, end, offset] = fmi.get_range(seed, 0);
 *     std::cout << "seed: " << seed << "\n";
 *     std::cout << "seed offset: " << offset << "\n";
 *     std::cout << "occurrence: " << end - beg << "\n";
 *     const auto offsets = fmi.get_offsets(beg, end);
 *     std::cout << "ref offsets: ";
 *     for (const auto offset : offsets) std::cout << offset << " ";
 *     std::cout << "\n";
 *   }
 * }
 * ```
 */
template<int SA_INTV = 1, typename size_type = std::uint32_t,
         SASorter Sorter = KISS1Sorter<size_type>>
class FMIndex {
 public:
  /**
   * Occ sampling interval, default value is 16.
   */
  const int OCC_INTV = 16;

  /**
   * The length of fixed suffix for lookup, default value is 14.
   */
  const int LOOKUP_LEN = 14;

  const int OCC1_INTV = 256;
  const int OCC2_INTV = OCC_INTV;

  using char_type = std::int8_t;

  /**
   * A compression vector which store the bwt.
   */
  DibitVector<std::uint8_t> bwt_;

  /**
   * A hierarchical sampled occurrence table.
   */
  std::pair<std::vector<std::array<size_type, 4>>,
            std::vector<std::array<std::uint8_t, 4>>>
    occ_;

  /**
   * A sampled suffix array.
   */
  Sorter::SA_t sa_;

  /**
   * A bit vector recording sampled suffix array.
   */
  detail::XbitVector<1, std::uint64_t, std::allocator<std::uint64_t>> b_;

  /**
   * prefix sum of b_ sampling interval, default value is 64.
   */
  const int B_OCC_INTV = 64;

  /**
   * A sampled prefix sum for b_.
   */
  std::vector<size_type> b_occ_;

  /**
   * A lookup table for fixed suffix query.
   */
  std::array<size_type, 4> cnt_{};
  size_type pri_{};
  std::vector<size_type> lookup_;

 protected:
  constexpr static auto cnt_table = [] {
    std::array<std::array<std::uint8_t, 4>, 256> cnt_table{};
    for (auto i = size_type{}; i < cnt_table.size(); i++)
      for (auto shift = size_type{}; shift < 8; shift += 2)
        cnt_table[i][i >> shift & 3u]++;
    return cnt_table;
  }();

  auto
  compute_occ(char_type c, size_type i) const {
    const auto occ1_beg = i / OCC1_INTV;
    const auto occ2_beg = i / OCC2_INTV;
    auto beg = occ2_beg * OCC2_INTV;
    auto cnt = size_type{};
    const auto pass_pri = c == 0 && beg <= pri_ && pri_ < i;
    const auto run = (i - beg) / 4;
    for (auto j = size_type{}; j < run; j++) {
      cnt += cnt_table[bwt_.data()[beg / 4]][c];
      beg += 4;
    }
    for (; beg < i; beg++)
      if (bwt_[beg] == c)
        cnt++;
    return occ_.first[occ1_beg][c] + occ_.second[occ2_beg][c] + cnt - pass_pri;
  }

  auto
  lf(char_type c, size_type i) const {
    return cnt_[c] + compute_occ(c, i);
  };

  auto
  compute_b_occ(size_type i) const {
    if constexpr (SA_INTV == 1)
      return i;
    else {
      const auto b_occ_beg = i / B_OCC_INTV;
      auto beg = b_occ_beg * B_OCC_INTV;
      auto cnt = size_type{};

      const auto run = (i - beg) / 64;
      for (auto j = size_type{}; j < run; j++) {
        cnt += std::popcount(b_.data()[beg / 64]);
        beg += 64;
      }

      auto mask = (1ull << i - beg) - 1;
      cnt += std::popcount(b_.data()[beg / 64] & mask);
      return b_occ_[b_occ_beg] + cnt;
    }
  }

  auto
  compute_sa(size_type i) const {
    if constexpr (SA_INTV == 1)
      return sa_[i];
    else {
      auto cnt = size_type{};
      while (not b_[i]) {
        i = lf(bwt_[i], i);
        cnt++;
      }
      return (sa_[compute_b_occ(i)] >> 2) * SA_INTV + cnt;
    }
  }

  auto
  compute_range(istring_view seed, size_type beg, size_type end,
                size_type stop_upper) const {
    while (!seed.empty()) {
      if (end - beg < stop_upper)
        break;
      beg = lf(seed.back(), beg);
      end = lf(seed.back(), end);
      seed.remove_suffix(1);
    }
    return std::array{beg, end, static_cast<size_type>(seed.size())};
  }

  auto
  build_lookup() {
    auto lookup_size = size_type{1} << LOOKUP_LEN * 2;
    lookup_.resize(lookup_size);
    lookup_.push_back(bwt_.size());

    const auto set_value = [this](const auto beg_i, const auto end_i) {
      const auto remain = (end_i - beg_i) % 2;
      const auto last_i = end_i - remain;
      for (auto i = beg_i; i < last_i; i += 2) {
        const auto seed = Codec::rhash(i, LOOKUP_LEN);
        const auto [beg, end, offset]
          = compute_range(seed, 0, lookup_.back(), 0);
        lookup_[i] = beg;
        lookup_[i + 1] = end;
      }
      if (remain) {
        const auto seed = Codec::rhash(last_i, LOOKUP_LEN);
        const auto [beg, end, offset]
          = compute_range(seed, 0, lookup_.back(), 0);
        lookup_[last_i] = beg;
      }
    };

    const auto thread_n = std::thread::hardware_concurrency();
    const auto batch_size = lookup_size / thread_n;
    auto workers = std::vector<std::thread>{};
    for (auto i = 0; i < thread_n; i++) {
      workers.emplace_back(set_value, i * batch_size, (i + 1) * batch_size);
    }
    set_value(thread_n * batch_size, lookup_size);
    for (auto& worker : workers) worker.join();
  }

  static auto
  validate_ref(istring_view s) {
    assert(std::all_of(std::execution::par_unseq, s.cbegin(), s.cend(),
                       [](auto c) { return c >= 0 && c <= 3; }));
  }

  void
  build_occ(istring_view ref, const auto& ori_sa) {
    auto& [occ1, occ2] = occ_;
    occ1.resize(ori_sa.size() / OCC1_INTV + 1);
    occ2.resize(ori_sa.size() / OCC2_INTV + 1);
#pragma omp parallel for
    for (auto beg = size_type{}; beg < ori_sa.size(); beg += OCC1_INTV) {
      auto& cnt = occ1[beg / OCC1_INTV];
      for (auto i = beg; i < beg + OCC1_INTV; i++) {
        if (i % OCC2_INTV == 0)
          for (int j = 0; j < 4; j++) occ2[i / OCC2_INTV][j] = cnt[j];
        if (i == ori_sa.size())
          break;
        const auto sa_v = ori_sa[i];
        if (sa_v != 0)
          cnt[ref[sa_v - 1]]++;
      }
    }

    for (auto& x : occ1) {
      for (int j = 0; j < 4; j++) {
        cnt_[j] += x[j];
        x[j] = cnt_[j] - x[j];
      }
    }

    auto sum = size_type{1};
    for (auto& x : cnt_) {
      sum += x;
      x = sum - x;
    }
  }

  void
  build_bwt(istring_view ref, const auto& ori_sa) {
    bwt_.resize(ori_sa.size());
    if constexpr (SA_INTV != 1)
      sa_.resize(ori_sa.size() / SA_INTV + 1);

#pragma omp parallel for
    for (auto beg = size_type{}; beg < ori_sa.size(); beg += 4) {
      const auto end = std::min((size_type)ori_sa.size(), beg + 4);
      for (auto i = beg; i < end; i++) {
        const auto sa_v = ori_sa[i];
        if (sa_v != 0) {
          bwt_[i] = ref[sa_v - 1];
        } else {
          bwt_[i] = 0;
          pri_ = i;
        }
      }
    }
  }

  void
  build_sa(istring_view ref, const auto& ori_sa) {
    if constexpr (SA_INTV == 1) {
      sa_ = ori_sa;
      return;
    }

    b_.resize(ori_sa.size());
    b_occ_.resize(ori_sa.size() / B_OCC_INTV + 1);
#pragma omp parallel for
    for (auto beg = size_type{}; beg < ori_sa.size(); beg += B_OCC_INTV) {
      auto& cnt = b_occ_[beg / B_OCC_INTV];
      for (auto i = beg; i < beg + B_OCC_INTV; i++) {
        if (i == ori_sa.size())
          break;
        const auto sa_v = ori_sa[i];
        b_[i] = sa_v % SA_INTV == 0;
        cnt += b_[i];
      }
    }

    auto sum = size_type{};
    for (auto& x : b_occ_) {
      sum += x;
      x = sum - x;
    }

    sa_.resize((ori_sa.size() + SA_INTV - 1) / SA_INTV);
#pragma omp parallel for
    for (auto beg = size_type{}; beg < ori_sa.size(); beg += B_OCC_INTV) {
      auto ptr = b_occ_[beg / B_OCC_INTV];
      for (auto i = beg; i < beg + B_OCC_INTV; i++) {
        if (i == ori_sa.size())
          break;
        const auto sa_v = ori_sa[i];
        if (sa_v % SA_INTV == 0) {
          sa_[ptr++] = sa_v;
        }
      }
    }
  }

 public:
  /**
   * Build index, initialize the core data structure from reference istring.
   *
   * Notice it will call Sorter::get_sa(...) to generate suffix array.
   */
  void
  build(istring_view ref) {
    SPDLOG_DEBUG("validate ref...");
    validate_ref(ref);

    const auto sort_len = 32u;
    SPDLOG_DEBUG("only build sa with prefix length: {}", sort_len);
    auto ori_sa = Sorter::get_suffix_array_dna(ref, sort_len);
    build(ref, ori_sa);
  }

  void
  build(istring_view ref, const auto& ori_sa) {
    // SPDLOG_DEBUG("validate sa...");
    // validate_sa(ref, ori_sa);

    SPDLOG_DEBUG("building FM-index begin...");
    SPDLOG_DEBUG("occ sampling interval: {}", OCC_INTV);
    SPDLOG_DEBUG("sa sampling interval: {}", SA_INTV);
    SPDLOG_DEBUG("lookup string length: {}", LOOKUP_LEN);
    if constexpr (SA_INTV != 1)
      SPDLOG_DEBUG("b occ inteval: {}", B_OCC_INTV);

    const auto thread_n = std::thread::hardware_concurrency();
    SPDLOG_DEBUG("using {} threads", thread_n);
    auto start = high_resolution_clock::now();
    auto end = high_resolution_clock::now();
    auto dur = duration_cast<seconds>(end - start);
    start = high_resolution_clock::now();

    build_occ(ref, ori_sa);
    build_bwt(ref, ori_sa);
    build_sa(ref, ori_sa);

    end = high_resolution_clock::now();
    dur = duration_cast<seconds>(end - start);
    SPDLOG_DEBUG("elapsed time: {} s.", dur.count());

    SPDLOG_DEBUG("computing {} suffix for for lookup...(using {} threads)",
                 (1ull << LOOKUP_LEN * 2), thread_n);
    start = high_resolution_clock::now();
    build_lookup();
    end = high_resolution_clock::now();
    dur = duration_cast<seconds>(end - start);
    SPDLOG_DEBUG("elapsed time: {} s.", dur.count());

    SPDLOG_DEBUG("validate lookup...");
    start = high_resolution_clock::now();

    assert(std::is_sorted(std::execution::par_unseq, lookup_.cbegin(),
                          lookup_.cend()));
    end = high_resolution_clock::now();
    dur = duration_cast<seconds>(end - start);
    SPDLOG_DEBUG("elapsed time: {} s.", dur.count());
  }

  auto
  get_offsets_traditional(size_type beg, size_type end) const {
    if constexpr (SA_INTV == 1)
      return std::span{&sa_[beg], end - beg};
    else {
      auto offsets = std::vector<size_type>{};
      offsets.reserve(end - beg);
      for (auto i = beg; i < end; i++) {
        offsets.push_back(compute_sa(i));
      }
      return offsets;
    }
  }

  /**
   * Compute the suffix array value according to the begin and end.
   * If sa_intv is 1, this can be done at O(1).
   */
  auto
  get_offsets(size_type beg, size_type end) const {
    if constexpr (SA_INTV == 1)
      return std::span{&sa_[beg], end - beg};
    else {
      auto offsets = std::vector<size_type>{};
      offsets.reserve(end - beg);

      constexpr auto MAX_DEPTH = SA_INTV;

      auto q = std::queue<std::tuple<size_type, size_type, int>>{};
      q.emplace(beg, end, 0);

      auto enqueue = [&q](auto beg, auto end, auto dep) {
        if (beg != end)
          q.emplace(beg, end, dep);
      };

      while (q.size() and offsets.size() < end - beg) {
        auto [cur_beg, cur_end, cur_dep] = q.front();
        q.pop();

        // add offset from sampled sa value
        const auto b_occ_cur_beg = compute_b_occ(cur_beg);
        const auto b_occ_cur_end = compute_b_occ(cur_end);
        for (auto i = b_occ_cur_beg; i < b_occ_cur_end; i++) {
          offsets.push_back(sa_[i] + cur_dep);
        }

        const auto nxt_dep = cur_dep + 1;
        if (nxt_dep == MAX_DEPTH)
          continue;

        if (cur_beg + 1 == cur_end) {
          const auto nxt_beg = lf(bwt_[cur_beg], cur_beg);
          const auto nxt_end = nxt_beg + 1;
          q.emplace(nxt_beg, nxt_end, nxt_dep);
        } else {
          [&]<auto... Idx>(std::index_sequence<Idx...>) {
            auto bg = std::array{lf(Idx, cur_beg)...};
            auto ed = std::array{lf(Idx, cur_end)...};
            (enqueue(bg[Idx], ed[Idx], nxt_dep), ...);
          }
          (std::make_index_sequence<4>{});
        }
      }
      return offsets;
    }
  }

  auto
  fmtree(istring_view seed) {
    auto [pre_beg, pre_end, pre_offs] = get_range(seed.substr(1), 0);
    auto [beg, end, offs] = get_range(seed.substr(0, 1), pre_beg, pre_end, 0);

    auto offsets = std::vector<size_type>{};
    offsets.reserve(end - beg);

    constexpr auto MAX_DEPTH = SA_INTV;

    auto q = std::queue<std::tuple<size_type, size_type, int>>{};
    q.emplace(beg, end, 0);

    auto enqueue = [&q](auto beg, auto end, auto dep) {
      if (beg != end)
        q.emplace(beg, end, dep);
    };

    while (q.size() and offsets.size() < end - beg) {
      auto [cur_beg, cur_end, cur_dep] = q.front();
      q.pop();

      // add offset from sampled sa value
      const auto b_occ_cur_beg = compute_b_occ(cur_beg);
      const auto b_occ_cur_end = compute_b_occ(cur_end);
      for (auto i = b_occ_cur_beg; i < b_occ_cur_end; i++) {
        offsets.push_back(sa_[i] + cur_dep);
      }

      const auto nxt_dep = cur_dep + 1;
      if (nxt_dep == MAX_DEPTH)
        continue;

      if (cur_beg + 1 == cur_end) {
        const auto nxt_beg = lf(bwt_[cur_beg], cur_beg);
        const auto nxt_end = nxt_beg + 1;
        q.emplace(nxt_beg, nxt_end, nxt_dep);
      } else {
        [&]<auto... Idx>(std::index_sequence<Idx...>) {
          auto bg = std::array{lf(Idx, cur_beg)...};
          auto ed = std::array{lf(Idx, cur_end)...};
          (enqueue(bg[Idx], ed[Idx], nxt_dep), ...);
        }
        (std::make_index_sequence<4>{});
      }
    }

    return offsets;
  }

  auto
  get_range(istring_view seed, size_type beg, size_type end,
            size_type stop_cnt = 0) const {
    if (end == beg || seed.empty())
      return std::array{beg, end, size_type{}};
    return compute_range(seed, beg, end, stop_cnt + 1);
  }

  /**
   * Get begin and end index of the uncompressed suffix array for the input
   * seed. Difference between begin and end is the occurrence count in
   * reference.
   *
   * @param seed Seed to search.
   * @param stop_cnt The remaining prefix length for the seed. Notice that when
   * the stop_cnt is -1, this value always be 0. The stop_cnt can be using to
   * early stop when occurrence count is not greater than the value.
   * Set to 0 to forbid early stop.
   * @return An array of begin, end index and the match stop position.
   */
  auto
  get_range(istring_view seed, size_type stop_cnt = 0) const {
    auto beg = size_type{};
    auto end = bwt_.size();
    if (seed.size() >= LOOKUP_LEN) {
      const auto key = Codec::hash(seed.substr(seed.size() - LOOKUP_LEN));
      beg = lookup_[key];
      end = lookup_[key + 1];
      seed.remove_suffix(LOOKUP_LEN);
    }
    return get_range(seed, beg, end, stop_cnt);
  }

  /**
   * Save index, utility for serialization.
   * The binary binary archive file size is same as the memory
   * occupation in run time.
   */
  auto
  save(std::ofstream& fout) const {
    const auto start = high_resolution_clock::now();
    fout.write(reinterpret_cast<const char*>(&cnt_), sizeof(cnt_));
    fout.write(reinterpret_cast<const char*>(&pri_), sizeof(pri_));
    SPDLOG_DEBUG("save bwt...");
    Serializer::save(fout, bwt_);
    SPDLOG_DEBUG("save occ...");
    Serializer::save(fout, occ_.first);
    Serializer::save(fout, occ_.second);
    SPDLOG_DEBUG("save sa...");
    Serializer::save(fout, sa_);
    SPDLOG_DEBUG("save lookup...");
    Serializer::save(fout, lookup_);

    if constexpr (SA_INTV != 1) {
      SPDLOG_DEBUG("save b_...");
      Serializer::save(fout, b_);
      SPDLOG_DEBUG("save b_occ_...");
      Serializer::save(fout, b_occ_);
    }
    const auto end = high_resolution_clock::now();
    const auto dur = duration_cast<seconds>(end - start);
    SPDLOG_DEBUG("elapsed time: {} s.", dur.count());
  }

  /**
   * Load index, utility for serialization.
   */
  auto
  load(std::ifstream& fin) {
    const auto start = high_resolution_clock::now();
    fin.read(reinterpret_cast<char*>(&cnt_), sizeof(cnt_));
    fin.read(reinterpret_cast<char*>(&pri_), sizeof(pri_));
    SPDLOG_DEBUG("load bwt...");
    Serializer::load(fin, bwt_);
    SPDLOG_DEBUG("load occ...");
    Serializer::load(fin, occ_.first);
    Serializer::load(fin, occ_.second);
    SPDLOG_DEBUG("load sa...");
    Serializer::load(fin, sa_);
    SPDLOG_DEBUG("load lookup...");
    Serializer::load(fin, lookup_);

    if constexpr (SA_INTV != 1) {
      SPDLOG_DEBUG("load b_...");
      Serializer::load(fin, b_);
      SPDLOG_DEBUG("load b_occ_...");
      Serializer::load(fin, b_occ_);
    }

    assert(fin.peek() == EOF);
    const auto end = high_resolution_clock::now();
    const auto dur = duration_cast<seconds>(end - start);
    SPDLOG_DEBUG("elapsed time: {} s.", dur.count());
  }

  bool
  operator==(const FMIndex& other) const = default;
};

}  // namespace biovoltron
