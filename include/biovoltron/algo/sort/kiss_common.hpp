#pragma once

#include <biovoltron/algo/sort/utils.hpp>

#include <algorithm>
#include <execution>
#include <span>

namespace biovoltron {

namespace kiss {

template <typename size_type>
auto get_S_bucket(
  const std::array<size_type, CHAR_SIZE> &count
) {
  auto bucket = count;
  std::inclusive_scan(
    bucket.begin(), bucket.end(),
    bucket.begin(),
    std::plus<>(), size_type{1}
  );
  return bucket;
}

template <typename size_type>
auto get_L_bucket(
  const std::array<size_type, CHAR_SIZE> &count
) {
  auto bucket = count;
  std::exclusive_scan(
    bucket.begin(), bucket.end(),
    bucket.begin(),
    size_type{1}
  );
  return bucket;
}

template <typename char_type, typename size_type>
void induced_S_block(
  const vector<char_type> &S,
  vector<size_type> &SA,
  std::array<size_type, CHAR_SIZE> &bucket,
  size_type block_beg,
  size_type block_end
) {
  while (block_end < block_beg) {
    auto sa_v = SA[block_beg];

    if (sa_v) {
      auto sa_u = sa_v - 1;
      auto cu = S[sa_u];
      auto cv = S[sa_v];
      if (cu <= cv)
        SA[--bucket[cu]] = sa_u;
    }

    block_beg--;
  }
}

template <typename char_type, typename size_type>
auto induced_S_block_prepare(
  const vector<char_type> &S,
  vector<size_type> &SA,
  ThreadState<size_type> &state
) {
  const size_t prefetch_distance = 32;
  auto &cache  = state.cache;
  auto &bucket = state.bucket;
  std::fill(state.bucket.begin(), state.bucket.begin() + CHAR_SIZE, 0);

  size_type i, j, count = 0;
  if (state.len == 0) return count;

  auto prefetchr_S = [&](size_type i) {
    auto s = SA[(ssize_t)i - prefetch_distance];
    const uint8_t *Ss = &S[s] - 1;
    prefetchr(s > 0 ? Ss : NULL);
    Ss--;
    prefetchr(s > 0 ? Ss : NULL);
  };

  auto prepare = [&](size_type i) {
    auto sa_v = SA[i];
    if (sa_v > 0) {
      auto sa_u = sa_v - 1;
      auto cv = S[sa_v];
      auto cu = S[sa_u];

      if (cu <= cv) {
        cache[count] = { cu, sa_u };
        bucket[cache[count++].symbol]++;
      }
    }
  };

  for (i = state.beg + state.len - 1, j = state.beg + prefetch_distance + 1; i >= j; i -= 2) {
    prefetchw(&SA[(ssize_t)i - 2 * prefetch_distance]);

    prefetchr_S(i - 0);
    prefetchr_S(i - 1);

    prepare(i - 0);
    prepare(i - 1);
  }

  for (j -= prefetch_distance + 1; i >= j; i--)
    prepare(i);

  return count;
}

template <typename char_type, typename size_type>
auto induced_S_block_update(
  const vector<char_type> &S,
  vector<size_type> &SA,
  ThreadState<size_type> &state
) {
  if (state.len == 0) return ;

  const size_t prefetch_distance = 32;
  auto &cache  = state.cache;
  auto &bucket = state.bucket;

  auto update = [&](size_type i) {
    SA[--bucket[cache[i].symbol]] = cache[i].index;
  };

  size_type i = 0, j = 0;
  if (state.m >= 3) {
    for (i = 0, j = state.m - 3; i < j; i += 4) {
      prefetchr(&cache[i + prefetch_distance]);

      update(i + 0);
      update(i + 1);
      update(i + 2);
      update(i + 3);
    }
  }

  for (j += std::min<size_type>(state.m, 3); i < j; i++)
    update(i);
}

template <typename char_type, typename size_type>
void induced_S_block(
  const vector<char_type> &S,
  vector<size_type> &SA,
  std::array<size_type, CHAR_SIZE> &bucket,
  size_type block_beg,
  size_type block_end,
  vector<ThreadState<size_type>> &states
) {
  auto num_threads = states.size();
  auto block_len = block_beg - block_end;

#pragma omp parallel num_threads(num_threads)
  {
    auto tid = omp_get_thread_num();
    auto &state = states[tid];
    distribute_workload<size_type>(block_len, state);
    state.beg += block_end + 1;

    // prepare
    state.m = induced_S_block_prepare(S, SA, state);

#pragma omp barrier

    // induced
#pragma omp master
    {
      for (size_t tid = num_threads - 1; ~tid; tid--) {
        auto &state = states[tid];
        for (auto c = 0; c < CHAR_SIZE; c++) {
          auto A = bucket[c];
          auto B = state.bucket[c];
          bucket[c] = A - B;
          state.bucket[c] = A;
        }
      }
    }

#pragma omp barrier

    // update
    induced_S_block_update(S, SA, state);
  }
}

template <typename char_type, typename size_type>
void induced_S(
  const vector<char_type> &S,
  vector<size_type> &SA,
  vector<ThreadState<size_type>> &states,
  const std::array<size_type, CHAR_SIZE> &count
) {
  auto n = S.size();
  auto num_threads = states.size();

  auto bucket = get_S_bucket(count);
  for (auto block_beg = n; block_beg > 0;) {
    if (SA[block_beg] == EMPTY<size_type>) {
      block_beg--;
      continue;
    }

    auto block_max_len = num_threads * (THREAD_CACHE_SIZE - 16 * num_threads);
    auto block_max_end = block_beg < block_max_len ? 0 : block_beg - block_max_len;
    auto block_end = block_beg - 1;
    while (block_end > block_max_end and SA[block_end] != EMPTY<size_type>)
      block_end--;

    auto block_len = block_beg - block_end;
    if (block_len < 32) {
      induced_S_block<char_type, size_type>(S, SA, bucket, block_beg, block_end);
    } else {
      induced_S_block<char_type, size_type>(S, SA, bucket, block_beg, block_end, states);
    }

    block_beg = block_end;
  }
}

template <typename char_type, typename size_type>
void induced_L_block(
  const vector<char_type> &S,
  vector<size_type> &SA,
  std::array<size_type, CHAR_SIZE> &bucket,
  size_type block_beg,
  size_type block_end
) {
  while (block_beg < block_end) {
    auto sa_v = SA[block_beg];

    if (sa_v) {
      auto sa_u = sa_v - 1;
      auto cv = sa_v == S.size() ? 0 : S[sa_v];
      auto cu = S[sa_u];
      SA[bucket[cu]++] = cu < cv ? EMPTY<size_type> : sa_u;
    }

    block_beg++;
  }
}

template <typename char_type, typename size_type>
auto induced_L_block_update(
  const vector<char_type> &S,
  vector<size_type> &SA,
  ThreadState<size_type> &state
) {
  if (state.len == 0) return ;

  const size_t prefetch_distance = 32;
  auto &cache  = state.cache;
  auto &bucket = state.bucket;

  auto update = [&](size_type i) {
    SA[bucket[cache[i].symbol]++] = cache[i].index;
  };

  size_type i, j;
  for (i = 0, j = state.m - 3; i < j; i += 4) {
    prefetchr(&cache[i + prefetch_distance]);

    update(i + 0);
    update(i + 1);
    update(i + 2);
    update(i + 3);
  }

  for (j += 3; i < j; i++)
    update(i);
}

template <typename char_type, typename size_type>
auto induced_L_block_prepare(
  const vector<char_type> &S,
  vector<size_type> &SA,
  ThreadState<size_type> &state
) {
  const size_t prefetch_distance = 32;
  auto &cache  = state.cache;
  auto &bucket = state.bucket;
  std::fill(state.bucket.begin(), state.bucket.begin() + CHAR_SIZE, 0);

  size_type i, j, count = 0;
  if (state.len == 0) return count;

  auto prefetchr_S = [&](size_type i) {
    auto s = SA[i + prefetch_distance];
    const uint8_t *Ss = &S[s] - 1;
    prefetchr(s > 0 ? Ss : NULL);
    Ss--;
    prefetchr(s > 0 ? Ss : NULL);
  };

  auto prepare = [&](size_type i) {
    auto sa_v = SA[i];
    if (sa_v != EMPTY<size_type> and sa_v > 0) {
      auto sa_u = sa_v - 1;
      auto cv = sa_v == S.size() ? 0 : S[sa_v];
      auto cu = S[sa_u];
      cache[count] = { cu, cu < cv ? EMPTY<size_type> : sa_u };
      bucket[cache[count++].symbol]++;
    }
  };

  for (i = state.beg, j = state.beg + state.len - prefetch_distance - 1; i < j; i += 2) {
    prefetchw(&SA[i + 2 * prefetch_distance]);

    prefetchr_S(i + 0);
    prefetchr_S(i + 1);

    prepare(i + 0);
    prepare(i + 1);
  }

  for (j += prefetch_distance + 1; i < j; i++)
    prepare(i);

  return count;
}

template <typename char_type, typename size_type>
void induced_L_block(
  const vector<char_type> &S,
  vector<size_type> &SA,
  std::array<size_type, CHAR_SIZE> &bucket,
  size_type block_beg,
  size_type block_end,
  vector<ThreadState<size_type>> &states
) {
  size_t num_threads = states.size();
  auto block_len = block_end - block_beg;

#pragma omp parallel num_threads(num_threads)
  {
    auto tid = omp_get_thread_num();
    auto &state = states[tid];
    distribute_workload<size_type>(block_len, state);
    state.beg += block_beg;

    // prepare
    state.m = induced_L_block_prepare(S, SA, state);

#pragma omp barrier

    // induced
#pragma omp master
    {
      for (auto tid = size_t{}; tid < num_threads; tid++) {
        auto &state = states[tid];
        for (auto c = 0; c < CHAR_SIZE; c++) {
          auto A = bucket[c];
          auto B = state.bucket[c];
          bucket[c] = A + B;
          state.bucket[c] = A;
        }
      }
    }

#pragma omp barrier

    // udpate
    induced_L_block_update(S, SA, state);
  }
}

template <typename char_type, typename size_type>
void induced_L(
  const vector<char_type> &S,
  vector<size_type> &SA,
  vector<ThreadState<size_type>> &states,
  const std::array<size_type, CHAR_SIZE> &count
) {
  auto n = S.size();
  auto num_threads = states.size();

  auto bucket = get_L_bucket(count);;
  for (auto block_beg = size_type{}; block_beg < n;) {
    if (SA[block_beg] == EMPTY<size_type>) {
      block_beg++;
      continue;
    }

    auto block_max_len = num_threads * (THREAD_CACHE_SIZE - 16 * num_threads);
    auto block_max_end = std::min<size_type>(n + 1, block_beg + block_max_len);
    auto block_end = block_beg + 1;
    while (block_end < block_max_end and SA[block_end] != EMPTY<size_type>)
      block_end++;

    auto block_len = block_end - block_beg;
    if (block_len < 32)
      induced_L_block<char_type, size_type>(S, SA, bucket, block_beg, block_end);
    else
      induced_L_block<char_type, size_type>(S, SA, bucket, block_beg, block_end, states);

    block_beg = block_end;
  }
}

template <typename size_type>
void induced_clear(
  vector<size_type> &SA,
  const std::array<size_type, CHAR_SIZE> &count,
  const std::array<size_type, CHAR_SIZE> &lms_count
) {
  auto bucket = get_S_bucket(count);
  for (auto i = 0; i < CHAR_SIZE; i++) {
    auto beg = bucket[i] - lms_count[i];
    auto end = bucket[i];
    if (beg == end)
      continue;

    std::fill(&SA[beg], &SA[end], EMPTY<size_type>);
  }
}

template <typename char_type, typename size_type>
void induced_sort(
  const vector<char_type> &S,
  vector<size_type> &SA,
  vector<ThreadState<size_type>> &states
) {
  auto num_threads = states.size();
  auto     count = std::array<size_type, CHAR_SIZE>{};
  auto lms_count = std::array<size_type, CHAR_SIZE>{};
  for (auto tid = size_t{}; tid < num_threads; tid++) {
    auto &bucket = states[tid].bucket;
    for (auto i = 0; i < CHAR_SIZE; i++) {
          count[i] += bucket[       4 * CHAR_SIZE + i];
      lms_count[i] += bucket[LMS_TYPE * CHAR_SIZE + i];
    }
  }

  induced_L(S, SA, states, count);
  induced_clear(SA, count, lms_count);
  induced_S(S, SA, states, count);
}


template <typename char_type, typename size_type>
void put_lms_suffix(
  const vector<char_type> &S,
  vector<size_type> &SA,
  const auto &SA1,
  vector<ThreadState<size_type>> &states
) {
  auto num_threads = states.size();
  auto     count = std::array<size_type, CHAR_SIZE>{};
  auto lms_count = std::array<size_type, CHAR_SIZE>{};
  for (auto tid = size_t{}; tid < num_threads; tid++) {
    auto &bucket = states[tid].bucket;
    for (auto i = 0; i < CHAR_SIZE; i++) {
          count[i] += bucket[       4 * CHAR_SIZE + i];
      lms_count[i] += bucket[LMS_TYPE * CHAR_SIZE + i];
    }
  }

  std::inclusive_scan(count.begin(), count.end(), count.begin(), std::plus<>(), size_type{1});

  size_type m = SA1.size();
  size_type j = SA.size();

  for (auto c = CHAR_SIZE - 1; ~c; c--) {
    auto num = lms_count[c];
    if (num == 0)
      continue;

    auto i = count[c];
    if (j - i > 0)
      std::fill(SA.begin() + i, SA.begin() + j, EMPTY<size_type>);

    std::memmove(&SA[j = (i - num)], &SA1[m -= num], sizeof(size_type) * num);
  }

  std::fill(SA.begin() + 1, SA.begin() + j, EMPTY<size_type>);
}

template <typename char_type, typename size_type, typename array_type>
auto get_lms(
  const vector<char_type> &S,
  array_type &SA,
  ThreadState<size_type> &state
) {
  // edge case
  if (state.len == 0)
    return size_type{};

  size_type n = S.size();
  auto end = state.beg + state.len, ptr = end;

  // determine the suffix type of S[end - 1]
  //  end-1   ptr
  //     v     v
  // -----|---------|
  //     ^     ^
  //    c0     c1
  auto c0 = S[end - 1], c1 = EMPTY<char_type>;
  while (ptr < n and c0 == (c1 = S[ptr])) ptr++;
  uint8_t type = ptr == n ? (bool)L_TYPE : c0 < c1;

  auto check_lms = [&](size_type i) {
    c1 = S[i - 1];

    // ------||
    //     ^^
    //    c1c0
    //
    // case1: c1 < c0, S_TYPE
    // case2: c1 > c0, L_TYPE
    // case3: c1 = c0, TYPE(c0)
    type = (type << 1) + (c1 == c0 ? (type & 1) : c1 < c0);

    auto curr_type = type & SUFFIX_TYPE_MASK;
    state.bucket[curr_type * CHAR_SIZE + c0]++;

    SA[ptr] = i;
    ptr += curr_type == LMS_TYPE;

    std::swap(c0, c1);
  };

  std::ranges::fill(state.bucket, 0);
  ptr = state.beg / 2;
  auto bucket = std::span{state.bucket.data() + 4 * CHAR_SIZE, CHAR_SIZE};
  for (auto i = end - 1; i >= state.beg + 1; i--) {
    bucket[c0]++;
    check_lms(i);
  }

  // no possible to be lms when i == 0
  bucket[c0]++;
  if (state.beg)
    check_lms(state.beg);

  return ptr - state.beg / 2;
}

template <typename char_type, typename size_type, typename array_type>
auto get_lms(
  const vector<char_type> &S,
  array_type &SA,
  vector<ThreadState<size_type>> &states
) {
  auto n = S.size();
  auto num_threads = states.size();

#pragma omp parallel num_threads(num_threads)
  {
    auto tid = omp_get_thread_num();
    auto &state = states[tid];
    distribute_workload<size_type>(n, state);
    state.m = get_lms(S, SA, state);
  }

  auto m = size_type{};
  for (auto tid = size_t{}; tid < num_threads; tid++) {
    auto &state = states[tid];

    // no need to move data when tid == 0 and no lms suffix
    if (tid and state.m)
      std::copy(
        SA.data() + state.beg / 2,
        SA.data() + state.beg / 2 + state.m,
        SA.data() + m
      );
    std::reverse(
      SA.data() + m,
      SA.data() + m + state.m
    );
    m += state.m;
  }

  return SA[m++] = n, m;
}

} // namespace biovoltron

} // namespace kiss
