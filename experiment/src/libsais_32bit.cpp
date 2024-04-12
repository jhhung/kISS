#define TBB_SUPPRESS_DEPRECATED_MESSAGES 1
#define SPDLOG_ACTIVE_LEVEL SPDLOG_LEVEL_DEBUG

#include "utils.hpp"
#include <libsais.h>
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

auto run(
  auto &S,
  size_t num_threads
) {
  auto sw = spdlog::stopwatch{};

  int32_t *SA = (int32_t *)calloc((S.size() + 1), sizeof(int32_t));
  const uint8_t *S_data = reinterpret_cast<const uint8_t *>(S.data());
  auto ret = libsais_omp(S_data, SA, S.size(), 0, NULL, num_threads);

  SPDLOG_DEBUG("libsais elapsed {}", sw);
  SPDLOG_DEBUG("Memory usage {}", getPeakRSS());
  free(SA);

  return ret;
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
  auto ret = run(S, num_threads);
  return ret;
}