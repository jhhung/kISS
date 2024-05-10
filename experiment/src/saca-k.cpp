#define TBB_SUPPRESS_DEPRECATED_MESSAGES 1
#define SPDLOG_ACTIVE_LEVEL SPDLOG_LEVEL_DEBUG

#include "utils.hpp"
#include <omp.h>
#include <spdlog/spdlog.h>
#include <spdlog/stopwatch.h>
#include <iostream>
#include <stdlib.h>
#include <sys/resource.h>

size_t getPeakRSS() {
  struct rusage rusage;
  getrusage(RUSAGE_SELF, &rusage);
  return (size_t)(rusage.ru_maxrss * 1024L);
}

void SACA_K(unsigned char *s, unsigned int *SA, unsigned int n, unsigned int K,
            unsigned int m, int level);

auto run(
  auto &S,
  size_t num_threads
) {
  // Move sequence to s_ch
  unsigned int n = S.size() + 1;
  unsigned char *s_ch = new unsigned char[n];

#pragma omp parallel for
  for (auto i = uint32_t{0}; i < S.size(); i++) {
    s_ch[i] = S[i];
  }
  S[n - 1] = 0;
  S.resize(0);
  S.shrink_to_fit();

  auto sw = spdlog::stopwatch{};

  unsigned int *SA = new unsigned int[n];
  SACA_K(s_ch, SA, n, 128, n, 0);

  SPDLOG_DEBUG("saca-k elapsed {}", sw);
  SPDLOG_DEBUG("Memory usage {}", getPeakRSS());
  
  delete[] SA;
  delete[] s_ch;
}

int main(int argc, char *argv[]) {
  spdlog::set_level(spdlog::level::debug);
  if (argc - 1 != 3) {
    std::cout << "Usage: " << argv[0] << " <fn.fa> <num_threads> <is_general>\n";
    exit(1);
  }

  auto S = std::string{};
  if (atoi(argv[3]) == 1) { // general sorting
    S = read_fasta(argv[1]);
  } else {
    S = read_fasta_dna(argv[1]);
  }

  size_t num_threads = atoi(argv[2]);
  run(S, num_threads);
}