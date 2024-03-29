#pragma once

#include "structs.hpp"
#include "constant.hpp"

#include <string>
#include <fstream>
#include <iostream>
#include <chrono>
#include <ranges>
#include <iomanip>
#include <algorithm>
#include <omp.h>

#define prefetchr(address) __builtin_prefetch((const void *)(address), 0, 3)
#define prefetchw(address) __builtin_prefetch((const void *)(address), 1, 3)

auto read_fasta(std::string filename) {
  auto fin = std::ifstream(filename);
  auto line = std::string{};
  auto seq = sais::vector<uint8_t>{};
  while (std::getline(fin, line)) {
    if (line[0] == '>')
      continue;
    for (auto &c : line) {
      seq.push_back(std::toupper(c));
    }
  }
  return seq;
}

namespace sais {

template <typename char_type, typename size_type>
void print(const vector<char_type> &S, const vector<size_type> &SA) {
  std::cout << std::setw(4) << "idx:";
  for (auto i = size_type{}; i < S.size() + 1; i++)
    std::cout << std::setw(4) << i;
  std::cout << '\n';

  std::cout << std::setw(4) << "S:";
  for (auto &v : S)
    std::cout << std::setw(4) << +v;
  std::cout << std::setw(4) << "$";
  std::cout << '\n';

  auto type = std::vector<char>{S_TYPE, L_TYPE};
  for (auto i = 1; i < S.size(); i++) {
    if (S[S.size() - i] == S[S.size() - i - 1])
      type.push_back(type.back());
    else
      type.push_back(S[S.size() - i - 1] < S[S.size() - i]);
  }
  std::reverse(type.begin(), type.end());

  std::cout << std::setw(4) << "T:";
  for (auto &v : type)
    std::cout << std::setw(4) << "LS"[v];
  std::cout << '\n';

  std::array<size_type, 128> wcnt{}, wlcnt{}, wscnt{};
  for (auto i = size_type{}; i < S.size(); i++) {
    auto c = S[i];
    wcnt[c]++;
    if (type[i])
      wscnt[c]++;
    else
      wlcnt[c]++;
  }

  std::cout << std::setw(4) << "BC:";
  std::cout << std::setw(4) << "$";
  for (auto c = 0; c < 128; c++)
    while (wcnt[c]--)
      std::cout << std::setw(4) << +c;
  std::cout << '\n';

  std::cout << std::setw(4) << "BT:";
  std::cout << std::setw(4) << "S";
  for (auto c = 0; c < 128; c++) {
    while (wlcnt[c]--)
      std::cout << std::setw(4) << "L";
    while (wscnt[c]--)
      std::cout << std::setw(4) << "S";
  }
  std::cout << '\n';

  std::cout << std::setw(4) << "SA:";
  for (auto &v : SA) {
    if (v == EMPTY<size_type>)
      std::cout << std::setw(4) << "-";
    else
      std::cout << std::setw(4) << +v;
  }
  std::cout << std::endl;
}

void check(const auto &S, const auto &SA, size_t sort_len) {
  size_t n = S.size();
  auto ok = true;
#pragma omp parallel for
  for (auto i = size_t{0}; i < n; i++) {
    auto ri = std::ranges::subrange(S.begin() + SA[i + 0], S.begin() + std::min(n, SA[i + 0] + sort_len));
    auto rj = std::ranges::subrange(S.begin() + SA[i + 1], S.begin() + std::min(n, SA[i + 1] + sort_len));
    if (SA[i + 0] > n or SA[i + 1] > n) {
      std::cout << "wrong value" << std::endl;
      ok = false;
    }
    if (std::ranges::lexicographical_compare(rj, ri)) {
      std::cout << "wrong order\n";
      std::cout << "SA[i + 0] = " << SA[i + 0] << " "; for (auto v : ri) std::cout << +v; std::cout << '\n';
      std::cout << "SA[i + 1] = " << SA[i + 1] << " "; for (auto v : rj) std::cout << +v; std::cout << '\n';
      std::cout << std::flush;
      ok = false;
    }
  }

  if (ok == false) {
    for (auto &v : S)
      std::cout << +v;
    std::cout << '\n';
    std::cout << std::flush;
    exit(0);
  }
}

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

}

