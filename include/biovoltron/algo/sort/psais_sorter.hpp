#pragma once

#include <biovoltron/utility/istring.hpp>
#include <biovoltron/algo/sort/psais_sorter.hpp>
#include <biovoltron/algo/sort/core/sorter.hpp>
#include <execution>
#include <iostream>

#include <thread>
#include <vector>
#include <string>
#include <array>
#include <limits>
#include <numeric>
#include <iomanip>
#include <future>
#include <ranges>
#include <execution>
#include <omp.h>

#include <boost/core/noinit_adaptor.hpp>
#include <biovoltron/utility/istring.hpp>
#include <biovoltron/utility/thread_pool.hpp>
#include <biovoltron/container/xbit_vector.hpp>

namespace biovoltron {

namespace psais::detail {

template<typename T>
std::reference_wrapper<std::remove_reference_t<T>> wrapper(T& t) { return std::ref(t); }

template<typename T>
T&& wrapper(std::remove_reference_t<T>&& t) { return std::move(t); }

} // psais::detail

namespace psais::utility {

template <typename Func, typename ... Args>
void parallel_do (
  std::integral auto n_jobs,
  std::integral auto n_threads,
  Func&& func,
  Args&& ...args
) {
  std::vector<std::jthread> threads;
  threads.reserve(n_jobs);
  auto counts = n_jobs / n_threads;
  auto remain = n_jobs % n_threads;

  for (decltype(n_jobs) tid = 0, L = 0; tid < n_threads; tid++) {
    auto block_size = counts + (tid < remain);
    if (block_size == 0) break;

    auto R = std::min(n_jobs, L + block_size);
    threads.emplace_back(std::forward<Func>(func), L, R, tid, psais::detail::wrapper<Args>(args)...);
    L = R;
  }
}

// #parallel_init
template <std::ranges::input_range R>
void parallel_init(R&& r, const auto& value) {
  std::fill(std::execution::par_unseq, r.begin(), r.end(), value);
}

// #parallel_prefix_sum
template <typename Vec>
void parallel_prefix_sum(
  Vec& vec,
  std::integral auto n_threads
) {
  auto n_jobs = vec.size();

  using T = Vec::value_type;

  std::vector<T, boost::noinit_adaptor<std::allocator<T>>> block_prefix_sum(n_jobs);
  std::vector<T> block_last_pos(n_threads, 0);
  biovoltron::psais::utility::parallel_do(vec.size(), (size_t) n_threads,
    [&](auto L, auto R, auto tid) {
      if (L == R)
        return ;

      block_prefix_sum[L] = vec[L];
      for (auto i = L + 1; i < R; i++)
        block_prefix_sum[i] = vec[i] + block_prefix_sum[i - 1];
      block_last_pos[tid] = R - 1;
    }
  );

  for (decltype(n_threads) i = 1; i < n_threads; i++) {
    if (block_last_pos[i] == 0)
      break;
    block_prefix_sum[block_last_pos[i]] += block_prefix_sum[block_last_pos[i - 1]];
  }

  biovoltron::psais::utility::parallel_do(n_jobs, (size_t) n_threads,
    [&](auto L, auto R, auto tid) {
      if (L == R)
        return ;

      auto offset = (tid == 0 ? 0 : block_prefix_sum[L - 1]);
      for (auto i = L; i < R - 1; i++)
        vec[i] = offset + block_prefix_sum[i];
      vec[R - 1] = block_prefix_sum[R - 1];
    }
  );
}

// #parallel_take_if
template <typename Compare, typename Project>
auto parallel_take_if(
  std::integral auto n_jobs,
  std::integral auto n_threads,
  std::ranges::random_access_range auto& ret,
  Compare &&compare,
  Project &&project
) {
  std::vector<size_t> buf(n_threads, 0);

  biovoltron::psais::utility::parallel_do(n_jobs, n_threads,
    [&](auto L, auto R, auto tid) {
      for (auto i = L; i < R; i++)
        buf[tid] += !!compare(i);
    }
  );

  std::vector<size_t> offset(n_threads + 1, 0);
  for (decltype(n_threads) i = 0; i < n_threads; i++)
    offset[i + 1] = offset[i] + buf[i];

  if constexpr (requires { ret.resize(0); }) {
    if (ret.size() < offset.back()) {
      ret.resize(offset.back());
    }
  }

  biovoltron::psais::utility::parallel_do(n_jobs, n_threads,
    [&](auto L, auto R, auto tid) {
      for (size_t i = L; i < R; i++)
        if (compare(i))
          ret[offset[tid]++] = project(i);
    }
  );
}

} //namespace psais::utility

namespace psais::detail {

#define L_TYPE 0
#define S_TYPE 1
#define NUM_THREADS std::thread::hardware_concurrency()
#define INDUCE_NUM_THREADS 16u

constexpr auto BLOCK_SIZE = 1u << 20;

template <typename T>
using NoInitVector = std::vector<T, boost::noinit_adaptor<std::allocator<T>>>;

using TypeVector = biovoltron::detail::XbitVector<1, std::uint8_t, std::allocator<std::uint8_t>>;

template <typename T>
constexpr auto EMPTY = std::numeric_limits<T>::max();

bool is_LMS(const TypeVector &T, auto i) {
  return i > 0 and T[i - 1] == L_TYPE and T[i] == S_TYPE;
}

// #name_substr
template<typename size_type>
auto name_substr(
  const std::ranges::random_access_range auto &S,
  const TypeVector &T,
  const NoInitVector<size_type> &SA,
  std::ranges::random_access_range auto S1,
  size_t sort_len
) {
  auto is_same_substr = [&S, &T] (auto x, auto y, auto k) {
    do {
      k--;
      if (S[x++] != S[y++]) return false;
    } while (!is_LMS(T, x) and !is_LMS(T, y) and k);

    return k and is_LMS(T, x) and is_LMS(T, y) and S[x] == S[y];
  };

  size_type n = (size_type)S.size();
  auto SA1 = NoInitVector<size_type>{};
  biovoltron::psais::utility::parallel_take_if(n, NUM_THREADS, SA1,
    [&](size_type i) { return is_LMS(T, SA[i]); },
    [&](size_type i) { return SA[i]; }
  );

  size_type n1 = (size_type)SA1.size();

  auto is_same = NoInitVector<size_type>(n1);
  biovoltron::psais::utility::parallel_init(is_same, 0);

  {
    auto result = std::vector<std::future<void>>{};
    result.reserve(n1 / BLOCK_SIZE + 1);

    auto pool = biovoltron::ThreadPool(NUM_THREADS);
    for (size_type x = 1; x < n1; x += BLOCK_SIZE) {
      size_type L = x, R = std::min(n1, L + BLOCK_SIZE);
      result.push_back(
        pool.enqueue(
          [&](size_type l, size_type r) {
            for (size_type i = l; i < r; i++)
              is_same[i] = not is_same_substr(SA1[i - 1], SA1[i], sort_len);
          }, L, R
        )
      );
    }

    for (auto &f : result) {
      f.get();
    }
  }

  biovoltron::psais::utility::parallel_prefix_sum(is_same, NUM_THREADS);

  NoInitVector<size_type> name(n);
  biovoltron::psais::utility::parallel_init(name, EMPTY<size_type>);

  biovoltron::psais::utility::parallel_do(n1, NUM_THREADS,
    [&](size_type L, size_type R, size_type) {
      for (size_type i = L; i < R; i++)
        name[SA1[i]] = is_same[i];
    }
  );

  biovoltron::psais::utility::parallel_take_if(n, NUM_THREADS, S1,
    [&](size_type i) { return name[i] != EMPTY<size_type>; },
    [&](size_type i) { return name[i]; }
  );

  return is_same.back();
}

// #induce_sort

// ##put_lms
template<typename size_type>
auto put_lms(
  const std::ranges::random_access_range auto &S,
  const NoInitVector<size_type> &LMS,
  const NoInitVector<size_type> &SA1,
  const NoInitVector<size_type> &BA,
  NoInitVector<size_type> &SA
) {
  size_type n1 = (size_type)SA1.size();
  size_type K = (size_type)BA.size() - 1;

  NoInitVector<size_type> S1(n1);
  std::transform(std::execution::par_unseq, SA1.begin(), SA1.end(), S1.begin(),
    [&LMS](auto &x) { return LMS[x]; }
  );

  NoInitVector<size_type> local_BA(1ull * (K + 1) * NUM_THREADS);
  biovoltron::psais::utility::parallel_do(NUM_THREADS, NUM_THREADS,
    [&](size_type, size_type, size_type tid) {
      size_type *ptr = local_BA.data() + 1ull * tid * (K + 1);
      for (size_type i = 0; i < K + 1; i++)
        ptr[i] = 0;
    }
  );

  biovoltron::psais::utility::parallel_do(n1, NUM_THREADS,
    [&](size_type L, size_type R, size_type tid) {
      size_type *ptr = local_BA.data() + 1ull * tid * (K + 1);
      for (size_type i = L; i < R; i++) {
        size_type idx = S1[i];
        ptr[S[idx]]++;
      }
    }
  );

  biovoltron::psais::utility::parallel_do(K + 1, NUM_THREADS,
    [&](size_type L, size_type R, size_type) {
      for (size_type i = NUM_THREADS - 2; ~i; i--) {
        auto *w_ptr = local_BA.data() + 1ull * (i    ) * (K + 1);
        auto *r_ptr = local_BA.data() + 1ull * (i + 1) * (K + 1);
        for (size_type j = L; j < R; j++)
          w_ptr[j] += r_ptr[j];
      }
    }
  );

  biovoltron::psais::utility::parallel_do(n1, NUM_THREADS,
    [&](size_type L, size_type R, size_type tid) {
      auto *ptr = local_BA.data() + 1ull * tid * (K + 1);
      for (size_type i = L; i < R; i++) {
        size_type idx = S1[i];
        size_type offset = ptr[S[idx]]--;
        SA[BA[S[idx]] - offset] = idx;
      }
    }
  );
}

// ##prepare
template<typename size_type, typename char_type>
void prepare(
  const size_t L,
  const std::ranges::random_access_range auto &S,
  NoInitVector<size_type> &SA,
  const TypeVector &T,
  NoInitVector<std::pair<char_type, uint8_t>> &RB
) {
  if (L >= SA.size()) return;
  decltype(L) R = std::min(SA.size(), L + BLOCK_SIZE);

  #pragma omp parallel for num_threads(INDUCE_NUM_THREADS / 2)
  for (auto i = L; i < R; i++) {
    auto induced_idx = SA[i] - 1;

    if (SA[i] == EMPTY<size_type> or SA[i] == 0) {
      RB[i - L] = {EMPTY<char_type>, 0};
    } else {
      RB[i - L] = {S[induced_idx], T[induced_idx]};
    }
  }
}

// ##update
template<typename size_type>
void update(
  const size_t L,
  const NoInitVector<std::pair<size_type, size_type>> &WB,
  NoInitVector<size_type> &SA
) {
  if (L >= SA.size()) return;
  decltype(L) R = std::min(SA.size(), L + BLOCK_SIZE);

  #pragma omp parallel for num_threads(INDUCE_NUM_THREADS / 2)
  for (auto i = L; i < R; i++) {
    auto& [idx, val] = WB[i - L];
    if (idx != EMPTY<size_type>) {
      SA[idx] = val;
    }
  }
}

// ##induce_impl
template<auto InduceType, typename size_type, typename char_type>
void induce_impl (
  const std::ranges::random_access_range auto &S,
  const TypeVector &T,
  const std::ranges::input_range auto &rng,
  const size_type L,
  NoInitVector<size_type> &SA,
  NoInitVector<std::pair<char_type, uint8_t>> &RB,
  NoInitVector<std::pair<size_type, size_type>> &WB,
  NoInitVector<size_type> &ptr
) {
  for (size_type i : rng) {
    auto induced_idx = SA[i] - 1;

    if (SA[i] != EMPTY<size_type> and SA[i] != 0) {
      auto chr = EMPTY<char_type>;
      if (auto [c, t] = RB[i - L]; c != EMPTY<char_type>) {
        if (t == InduceType) chr = c;
      } else if (T[induced_idx] == InduceType) {
        chr = S[induced_idx];
      }

      if (chr == EMPTY<char_type>) continue;

      bool is_adjacent;
      auto pos = ptr[chr];
      if constexpr (InduceType == L_TYPE) {
        ptr[chr] += 1;
        is_adjacent = pos < L + (BLOCK_SIZE << 1);
      } else {
        ptr[chr] -= 1;
        is_adjacent = pos + BLOCK_SIZE >= L;
      }

      // if pos is in adjacent block -> directly write it
      // otherwise, write it to WB
      if (is_adjacent) {
        SA[pos] = induced_idx;
        WB[i - L].first = EMPTY<size_type>;
      } else {
        WB[i - L] = {pos, induced_idx};
      }
    }
  }
}

// ##induce
template<auto InduceType, typename size_type, typename char_type>
void induce (
  const std::ranges::random_access_range auto &S,
  const TypeVector &T,
  NoInitVector<size_type> &SA,
  NoInitVector<std::pair<char_type, uint8_t>> &RBP,
  NoInitVector<std::pair<char_type, uint8_t>> &RBI,
  NoInitVector<std::pair<size_type, size_type>> &WBU,
  NoInitVector<std::pair<size_type, size_type>> &WBI,
  NoInitVector<size_type> &ptr
) {
  // views
  constexpr auto iter_view = [] {
    if constexpr (InduceType == L_TYPE) {
      return std::views::all;
    } else {
      return std::views::reverse;
    }
  }();

  size_type size = SA.size();
  auto blocks = std::views::iota(size_type(0), size)
    | std::views::filter([](size_type n) { return n % BLOCK_SIZE == 0; });

  // prepare for first block
  if constexpr (InduceType == L_TYPE) {
    prepare(0, S, SA, T, RBP);
  } else {
    prepare(size / BLOCK_SIZE * BLOCK_SIZE, S, SA, T, RBP);
  }

  auto pool = biovoltron::ThreadPool(2);
  auto stage = std::array<std::future<void>, 2>{};

  // pipeline
  for (size_type L : blocks | iter_view) {
    for (auto &s : stage) if (s.valid()) s.wait();
    RBI.swap(RBP);
    WBI.swap(WBU);

    // prepare && update
    size_type P_L = L + BLOCK_SIZE;
    size_type U_L = L - BLOCK_SIZE;
    if constexpr (InduceType == S_TYPE) {
      std::swap(P_L, U_L);
    }

    stage[0] = pool.enqueue(prepare<size_type, char_type, decltype(S)>, P_L,
        std::ref(S), std::ref(SA), std::ref(T), std::ref(RBP));
    stage[1] = pool.enqueue(update<size_type>, U_L,
        std::ref(WBU), std::ref(SA));

    // induce
    auto rng = std::views::iota(L, std::min(L + BLOCK_SIZE, size)) | iter_view;
    induce_impl<InduceType>(S, T, rng, L, SA, RBI, WBI, ptr);
  }
}

// ##induce_sort
template<typename size_type>
void induce_sort(
  const std::ranges::random_access_range auto &S,
  const TypeVector &T,
  const NoInitVector<size_type> &SA1,
  const NoInitVector<size_type> &LMS,
  NoInitVector<size_type> &BA,
  NoInitVector<size_type> &SA
) {
  using char_type = decltype(S.begin())::value_type;

  // induce LMS
  put_lms(S, LMS, SA1, BA, SA);

  // declare ptr, RBP, RBI, WBI, WBU
  NoInitVector<size_type> ptr(BA.size());
  NoInitVector<std::pair<char_type, uint8_t>> RBP(BLOCK_SIZE), RBI(BLOCK_SIZE);
  NoInitVector<std::pair<size_type, size_type>> WBU(BLOCK_SIZE), WBI(BLOCK_SIZE);

  // init buffer
  biovoltron::psais::utility::parallel_init(RBP, std::pair{EMPTY<char_type>, uint8_t(0)});
  biovoltron::psais::utility::parallel_init(RBI, std::pair{EMPTY<char_type>, uint8_t(0)});
  biovoltron::psais::utility::parallel_init(WBU, std::pair{EMPTY<size_type>, EMPTY<size_type>});
  biovoltron::psais::utility::parallel_init(WBI, std::pair{EMPTY<size_type>, EMPTY<size_type>});

  // induce L
  ptr[0] = 0;
  std::transform(std::execution::par_unseq, BA.begin(), BA.end() - 1, ptr.begin() + 1,
    [](size_type &b) { return b; }
  );
  induce<L_TYPE>(S, T, SA, RBP, RBI, WBU, WBI, ptr);

  // init buffer
  biovoltron::psais::utility::parallel_init(RBP, std::pair{EMPTY<char_type>, uint8_t(0)});
  biovoltron::psais::utility::parallel_init(RBI, std::pair{EMPTY<char_type>, uint8_t(0)});
  biovoltron::psais::utility::parallel_init(WBU, std::pair{EMPTY<size_type>, EMPTY<size_type>});
  biovoltron::psais::utility::parallel_init(WBI, std::pair{EMPTY<size_type>, EMPTY<size_type>});

  // clean S_TYPE
  std::for_each(std::execution::par_unseq, SA.begin() + 1, SA.end(),
    [&T](size_type &idx) {
      if (idx != EMPTY<size_type> and T[idx] == S_TYPE) {
        idx = EMPTY<size_type>;
      }
    }
  );

  // induce S
  std::transform(std::execution::par_unseq, BA.begin(), BA.end(), ptr.begin(),
    [](auto &bucket) { return bucket - 1; }
  );
  induce<S_TYPE>(S, T, SA, RBP, RBI, WBU, WBI, ptr);
}

// #preprocess

// ##get_type
template<typename size_type>
auto get_type(const std::ranges::random_access_range auto &S) {
  auto T = TypeVector(S.size());
  std::vector<size_type> same_char_suffix_len(NUM_THREADS, 0);
  std::vector<size_type> block_size(NUM_THREADS, 0);
  std::vector<size_type> block_left(NUM_THREADS, 0);

  auto cal_type = [&](auto x) -> bool {
    if (x + 1 == S.size())
      return S_TYPE;
    auto x1 = S[x], x2 = S[x + 1];
    if (x1 < x2)
      return S_TYPE;
    else if (x1 > x2)
      return L_TYPE;
    else
      return T[x + 1];
  };

  T[S.size() - 1] = S_TYPE;
  biovoltron::psais::utility::parallel_do((S.size() + 7) / 8, NUM_THREADS,
    [&](size_type L, size_type R, size_type tid) {
      if (L == R)
        return ;

      L = L * 8, R = std::min((size_type)S.size(), R * 8);
      if (R != S.size())
        T[R - 1] = cal_type(R - 1);

      same_char_suffix_len[tid] = 1;
      bool same = true;
      for (size_type i = R - L - 2; ~i; i--) {
        size_type x = L + i;
        T[x] = cal_type(x);

        if (S[x] != S[x + 1])
          same = false;

        if (same)
          same_char_suffix_len[tid]++;
      }

      block_size[tid] = R - L;
      block_left[tid] = L;
    }
  );

  std::vector<uint8_t> flip(NUM_THREADS, false);
  for (size_type i = NUM_THREADS - 2; ~i; i--) {
    if (block_size[i + 1] == 0)
      continue;

    size_type x1 = block_left[i + 1] - 1;
    size_type x2 = block_left[i + 1];
    // ...-|----|----|-...
    //        x1 x2

    if (S[x1] != S[x2])
      continue;

    uint8_t prev_left_type = T[x2];
    if (same_char_suffix_len[i + 1] == block_size[i + 1] and flip[i + 1])
      prev_left_type ^= 1;

    if (T[x1] != prev_left_type)
      flip[i] = true;
  }

  biovoltron::psais::utility::parallel_do((S.size() + 7) / 8, NUM_THREADS,
    [&](size_type L, size_type R, size_type tid) {
      if (not flip[tid])
        return ;

      L = L * 8, R = std::min((size_type)S.size(), R * 8);
      T[R - 1] = !T[R - 1];
      for (size_type i = R - L - 2; ~i; i--) {
        size_type x = L + i;
        if (S[x] != S[x + 1])
          return ;
        T[x] = !T[x];
      }
    }
  );

  return T;
}

// ##get_bucket
template<typename size_type>
auto get_bucket(const std::ranges::random_access_range auto &S, size_type K) {
  NoInitVector<size_type> local_BA(1ull * (K + 1) * NUM_THREADS);
  biovoltron::psais::utility::parallel_init(local_BA, 0);

  size_type n = S.size();
  biovoltron::psais::utility::parallel_do(n, NUM_THREADS,
    [&](size_type L, size_type R, size_type tid) {
      size_type *ptr = local_BA.data() + 1ull * tid * (K + 1);
      for (size_type i = L; i < R; i++)
        ptr[S[i]]++;
    }
  );

  auto BA = NoInitVector<size_type>(K + 1);
  biovoltron::psais::utility::parallel_init(BA, 0);

  biovoltron::psais::utility::parallel_do(K + 1, NUM_THREADS,
    [&](size_type L, size_type R, size_type) {
      for (size_type i = 0; i < NUM_THREADS; i++) {
        size_type *ptr = local_BA.data() + 1ull * i * (K + 1);
        for (size_type j = L; j < R; j++)
          BA[j] += ptr[j];
      }
    }
  );

  biovoltron::psais::utility::parallel_prefix_sum(BA, NUM_THREADS);

  return BA;
}

// ##get_lms
template<typename size_type>
auto get_lms(const TypeVector &T, const auto size) {
  auto LMS = NoInitVector<size_type>{};
  biovoltron::psais::utility::parallel_take_if(size, NUM_THREADS, LMS,
    [&](size_type i) { return is_LMS(T, i); },
    [ ](size_type i) { return i; }
  );
  return LMS;
}

// #suffix_array
template<typename size_type>
NoInitVector<size_type> suffix_array(
  const std::ranges::random_access_range auto &S,
  size_type K,
  size_t sort_len
) {
  // 1. get type && bucket array
  auto T = get_type<size_type>(S);
  auto BA = get_bucket(S, K);

  // 2. induce LMS-substring
  auto LMS = get_lms<size_type>(T, S.size());

  auto SA = NoInitVector<size_type>(S.size());
  biovoltron::psais::utility::parallel_init(SA, EMPTY<size_type>);

  auto SA1 = NoInitVector<size_type>(LMS.size());

  // iota SA1
  auto iota = std::views::iota(size_type(0), static_cast<size_type>(SA1.size()));
  std::transform(std::execution::par_unseq, iota.begin(), iota.end(), SA1.begin(),
    [](auto &idx) { return idx; }
  );

  induce_sort(S, T, SA1, LMS, BA, SA);

  auto S1 = std::ranges::subrange(SA.begin() + SA.size() - SA1.size(), SA.end());
  auto K1 = name_substr(S, T, SA, S1, sort_len);

  // 3. recursively solve LMS-suffix
  if (K1 + 1 == LMS.size()) {
    for (size_t i = 0; i < LMS.size(); i++) {
      SA1[S1[i]] = i;
    }
  } else {
    SA1 = suffix_array(S1, K1, sort_len >> 1);
  }

  // 4. induce orig SA
  biovoltron::psais::utility::parallel_init(SA, EMPTY<size_type>);
  induce_sort(S, T, SA1, LMS, BA, SA);

  return SA;
}

#undef L_TYPE
#undef S_TYPE
#undef NUM_THREADS
#undef INDUCE_NUM_THREADS

} // namespace psais::detail

namespace psais {

template <typename size_type>
auto suffix_array(istring_view s, size_t sort_len = std::string::npos) {
  auto ref = psais::detail::NoInitVector<uint8_t>(s.size() + 1);
  std::transform(std::execution::par_unseq, s.begin(), s.end(), ref.begin(),
    [](auto &c) { return c + 1; }
  );
  ref[s.size()] = 0;

#ifdef DEBUG
  auto bg = std::chrono::high_resolution_clock::now();
#endif
  auto sa = psais::detail::suffix_array(ref, 5u, sort_len);
#ifdef DEBUG
  auto ed = std::chrono::high_resolution_clock::now();
  std::cout << "psais " << (ed - bg).count() / 1e9 << "s. \n";
#endif
  return sa;
}

} // namespace psais

template<typename size_type = std::uint32_t>
struct PsaisSorter {
  static auto
  get_sa(istring_view ref, std::size_t sort_len = istring_view::npos) {
    auto sa = biovoltron::psais::suffix_array<size_type>(ref, sort_len);
    return std::vector<size_type>(sa.begin(), sa.end());
  }
};

}  // namespace biovoltron

