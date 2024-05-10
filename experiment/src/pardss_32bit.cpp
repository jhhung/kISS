#define TBB_SUPPRESS_DEPRECATED_MESSAGES 1
#define SPDLOG_ACTIVE_LEVEL SPDLOG_LEVEL_DEBUG

#include "utils.hpp"
#include <divsufsort.h>
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

template <typename num_type>
auto run(
  auto &S,
  size_t num_threads
) {
  auto sw = spdlog::stopwatch{};

  num_type size = S.size();
  num_type *SA = new num_type[size];
  divsufsort((sauchar_t *)S.data(), SA, size);

  SPDLOG_DEBUG("pardss elapsed {}", sw);
  SPDLOG_DEBUG("Memory usage {}", getPeakRSS());

  delete[] SA;
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
  run<int32_t>(S, num_threads);
}
