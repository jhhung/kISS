#pragma once

#include <biovoltron/algo/sort/structs.hpp>
#include <omp.h>

#define prefetchr(address) __builtin_prefetch((const void *)(address), 0, 3)
#define prefetchw(address) __builtin_prefetch((const void *)(address), 1, 3)

namespace biovoltron {

namespace kiss {

template <typename size_type>
struct alignas(64) ThreadState;

template <typename size_type>
void distribute_workload(size_type n, ThreadState<size_type> &state) {
  size_type thread_num  = omp_get_thread_num();
  size_type num_threads = omp_get_num_threads();
  size_type workload_len = (n / num_threads) & (-16);
  size_type beg = workload_len * thread_num;
  size_type len = thread_num + 1 == num_threads ? n - beg : workload_len;

  state.beg = beg;
  state.len = len;
}

template <typename size_type>
void distribute_workload_packed_DNA(size_type n, ThreadState<size_type> &state) {
  size_type thread_num  = omp_get_thread_num();
  size_type num_threads = omp_get_num_threads();
  size_type workload_len = (n / num_threads) & (-64);
  size_type beg = workload_len * thread_num;
  size_type len = thread_num + 1 == num_threads ? n - beg : workload_len;

  state.beg = beg;
  state.len = len;
}

} // namespace kiss

} // namespace biovoltron
