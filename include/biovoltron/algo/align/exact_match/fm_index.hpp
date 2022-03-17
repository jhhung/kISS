#pragma once

#include <biovoltron/algo/sort/core/sorter.hpp>
#include <biovoltron/algo/sort/psais_sorter.hpp>
#include <biovoltron/container/xbit_vector.hpp>
#include <biovoltron/utility/archive/serializer.hpp>
#include <biovoltron/utility/istring.hpp>
#include <chrono>
#include <span>
#include <thread>
#ifdef DEBUG
#include <iostream>
#endif

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
 * - hierarchical occurrence table: L1 occ occupy fixed `3.1Gb * 16 / 256 = 0.194Gb`
 *   plus L2 occ occupy `3.1Gb * 4 / occ_intv(16) = 0.775 Gb`.
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

template<
  typename size_type = std::uint32_t,
  SASorter Sorter = PsaisSorter<size_type>
>
class FMIndex {
 public:
  /**
   * Sa sampling interval, default value is 1.
   */
  static constexpr int SA_INTV = 1;

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
  std::vector<size_type> sa_;

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
  compute_sa(size_type i) const {
    auto cnt = size_type{};
    while (i % SA_INTV && i != pri_) {
      i = lf(bwt_[i], i);
      cnt++;
    }
    return i != pri_ ? sa_[i / SA_INTV] + cnt : cnt;
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
  compute_lookup() {
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
  validate(istring_view s) {
    assert(std::all_of(std::execution::par_unseq, s.cbegin(), s.cend(),
                       [](auto c) { return c >= 0 && c <= 3; }));
  }

 public:
  /**
   * Build index, initialize the core data structure from reference istring.
   *
   * Notice it will call Sorter::get_sa(...) to generate suffix array.
   */
  void 
  build(istring_view ref) {
#ifdef DEBUG
    std::cout << "build FM-index begin...\n";
    std::cout << "occ sampling interval: " << OCC_INTV << "\n";
    std::cout << "sa sampling interval: " << SA_INTV << "\n";
    std::cout << "lookup string length: " << LOOKUP_LEN << "\n";
#endif

#ifdef DEBUG
    std::cout << "validate ref...\n";
    auto start = high_resolution_clock::now();
#endif
    validate(ref);
#ifdef DEBUG
    auto end = high_resolution_clock::now();
    auto dur = duration_cast<seconds>(end - start);
    std::cout << "elapsed time: " << dur.count() << " s.\n";
#endif

    const auto sort_len = std::same_as<Sorter, PsaisSorter<size_type>> ? istring::npos : 256u;
#ifdef DEBUG
    std::cout << "only build sa with prefix length: " << sort_len << "\n";
    const auto thread_n = std::thread::hardware_concurrency();
    std::cout << "sa sort start...(using " << thread_n << " threads)\n";
    start = high_resolution_clock::now();
#endif
    auto ori_sa = Sorter::get_sa(ref, sort_len);
#ifdef DEBUG
    end = high_resolution_clock::now();
    dur = duration_cast<seconds>(end - start);
    std::cout << "elapsed time: " << dur.count() << " s.\n";
#endif

#ifdef DEBUG
    std::cout << "build bwt and sample occ...\n";
    start = high_resolution_clock::now();
#endif
    auto& [occ1, occ2] = occ_;
    occ1.resize(ori_sa.size() / OCC1_INTV + 2);
    occ2.resize(ori_sa.size() / OCC2_INTV + 1);

#pragma omp parallel for
    for (auto block_idx = size_type{}; block_idx < ori_sa.size(); block_idx += OCC1_INTV) {
      auto cnt = std::array<size_type, 4>{};
      for (auto offset = size_type{}; offset < OCC1_INTV; offset++) {
        auto i = block_idx + offset;
        if (i >= ori_sa.size())
          break;

        if (i % OCC2_INTV == 0)
          for (int j = 0; j < 4; j++)
            occ2[i / OCC2_INTV][j] = cnt[j];

        const auto sa_v = ori_sa[i];
        if (sa_v != 0) {
          const auto c = ref[sa_v - 1];
          cnt[c]++;
        }
      }
      occ1[block_idx / OCC1_INTV + 1] = cnt;
    }

    for (auto i = size_type{}; i < occ1.size(); i++) {
      for (int j = 0; j < 4; j++)
        cnt_[j] += occ1[i][j];
      occ1[i] = cnt_;
    }

    auto sum = size_type{1};
    for (auto& x : cnt_) {
      sum += x;
      x = sum - x;
    }

    bwt_.resize(ori_sa.size());
    sa_.resize(ori_sa.size() / SA_INTV + 1);
#pragma omp parallel for
    for (auto block_size = size_type{}; block_size < ori_sa.size(); block_size += 4) {
      for (int offset = 0; offset < 4; offset++) {
        auto i = block_size + offset;
        if (i >= ori_sa.size())
          break;
        const auto sa_v = ori_sa[i];
        if (sa_v != 0) {
          bwt_[i] = ref[sa_v - 1];
        } else {
          bwt_[i] = 0;
          pri_ = i;
        }
        if constexpr (SA_INTV != 1) {
          if (i % SA_INTV == 0)
            sa_[i / SA_INTV] = sa_v;
        }
      }
    }

    if constexpr (SA_INTV == 1)
      sa_.swap(ori_sa);
#ifdef DEBUG
    end = high_resolution_clock::now();
    dur = duration_cast<seconds>(end - start);
    std::cout << "elapsed time: " << dur.count() << " s.\n";
#endif

#ifdef DEBUG
    std::cout << "computing " << (1ull << LOOKUP_LEN * 2)
              << " suffix for for lookup...(using " << thread_n
              << " threads)\n";
    start = high_resolution_clock::now();
#endif
    compute_lookup();
#ifdef DEBUG
    end = high_resolution_clock::now();
    dur = duration_cast<seconds>(end - start);
    std::cout << "elapsed time: " << dur.count() << " s.\n";
#endif

#ifdef DEBUG
    std::cout << "validate lookup...\n";
    start = high_resolution_clock::now();
#endif

    assert(std::is_sorted(std::execution::par_unseq, lookup_.cbegin(),
                          lookup_.cend()));
#ifdef DEBUG
    end = high_resolution_clock::now();
    dur = duration_cast<seconds>(end - start);
    std::cout << "elapsed time: " << dur.count() << " s.\n";
#endif
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
      for (auto i = beg; i < end; i++) offsets.push_back(compute_sa(i));
      return offsets;
    }
  }

  auto
  get_range(istring_view seed, size_type beg, size_type end,
            size_type stop_cnt) const {
    if (end == beg || seed.empty())
      return std::array{beg, end, size_type{}};
    return compute_range(seed, beg, end, stop_cnt + 1);
  }

  /**
   * Get begin and end index of the uncompressed suffix array for the input seed.
   * Difference between begin and end is the occurrence count in reference.
   *
   * @param seed Seed to search.
   * @param stop_cnt The remaining prefix length for the seed. Notice that when 
   * the stop_cnt is -1, this value always be 0. The stop_cnt can be using to 
   * early stop when occurrence count is not greater than the value.
   * Set to -1 to forbid early stop.
   * @return An array of begin, end index and the match stop position.
   */
  auto
  get_range(istring_view seed, size_type stop_cnt) const {
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
#ifdef DEBUG
    const auto start = high_resolution_clock::now();
#endif
    fout.write(reinterpret_cast<const char*>(&cnt_), sizeof(cnt_));
    fout.write(reinterpret_cast<const char*>(&pri_), sizeof(pri_));
#ifdef DEBUG
    std::cout << "save bwt...\n";
#endif
    Serializer::save(fout, bwt_);
#ifdef DEBUG
    std::cout << "save occ...\n";
#endif
    Serializer::save(fout, occ_.first);
    Serializer::save(fout, occ_.second);
#ifdef DEBUG
    std::cout << "save sa...\n";
#endif
    Serializer::save(fout, sa_);
#ifdef DEBUG
    std::cout << "save lockup...\n";
#endif
    Serializer::save(fout, lookup_);
#ifdef DEBUG
    const auto end = high_resolution_clock::now();
    const auto dur = duration_cast<seconds>(end - start);
    std::cout << "elapsed time: " << dur.count() << " s.\n";
#endif
  }

  /**
   * Load index, utility for serialization.
   */
  auto
  load(std::ifstream& fin) {
#ifdef DEBUG
    const auto start = high_resolution_clock::now();
#endif
    fin.read(reinterpret_cast<char*>(&cnt_), sizeof(cnt_));
    fin.read(reinterpret_cast<char*>(&pri_), sizeof(pri_));
#ifdef DEBUG
    std::cout << "load bwt...\n";
#endif
    Serializer::load(fin, bwt_);
#ifdef DEBUG
    std::cout << "load occ...\n";
#endif
    Serializer::load(fin, occ_.first);
    Serializer::load(fin, occ_.second);
#ifdef DEBUG
    std::cout << "load sa...\n";
#endif
    Serializer::load(fin, sa_);
#ifdef DEBUG
    std::cout << "load lookup...\n";
#endif
    Serializer::load(fin, lookup_);
    assert(fin.peek() == EOF);
#ifdef DEBUG
    const auto end = high_resolution_clock::now();
    const auto dur = duration_cast<seconds>(end - start);
    std::cout << "elapsed time: " << dur.count() << " s.\n";
#endif
  }

  bool
  operator==(const FMIndex& other) const = default;
};

}  // namespace biovoltron
