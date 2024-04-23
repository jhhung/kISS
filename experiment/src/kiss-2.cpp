#define TBB_SUPPRESS_DEPRECATED_MESSAGES 1
#define SPDLOG_ACTIVE_LEVEL SPDLOG_LEVEL_DEBUG

#include <biovoltron/algo/sort/kiss2_core.hpp>
#include <biovoltron/algo/sort/utils.hpp>
#include <tbb/global_control.h>
#include "utils.hpp"
#include "omp.h"
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

template <typename size_type>
auto run(
  auto &S,
  size_t num_threads,
  int is_general,
  size_type k
) {
  auto S_input = biovoltron::kiss::vector<uint8_t>{S.begin(), S.end()};
  S.resize(0);
  S.shrink_to_fit();
  auto SA = biovoltron::kiss::vector<size_type>{};
  
  if (is_general) {
    auto sw = spdlog::stopwatch{};

    biovoltron::kiss::kiss2_suffix_array(S_input, SA, k, num_threads);

    SPDLOG_DEBUG("kiss-2 elapsed {}", sw);
    SPDLOG_DEBUG("Memory usage {}", getPeakRSS());

  } else {
    auto sw = spdlog::stopwatch{};

    biovoltron::kiss::kiss2_suffix_array_dna(S_input, SA, k, num_threads);

    SPDLOG_DEBUG("kiss-2 elapsed {}", sw);
    SPDLOG_DEBUG("Memory usage {}", getPeakRSS());
    
  }
}

int main(int argc, char *argv[]) {
  spdlog::set_level(spdlog::level::debug);
  if (argc - 1 != 4) {
    std::cout << "Usage: " << argv[0] << " <fn.fa> <num_threads> <is_general> <k>\n";
    exit(1);
  }

  auto S = std::string{};
  if (atoi(argv[3]) == 1) { // general sorting
    S = read_fasta(argv[1]);
  } else {
    S = read_fasta_dna(argv[1]);
#pragma omp parallel for
    for (auto i = uint32_t{0}; i < S.size(); i++) {
      S[i]--;
    }
  }

  size_t num_threads = atoi(argv[2]);
  auto global_control = tbb::global_control(
    tbb::global_control::max_allowed_parallelism,
    num_threads
  );
  run<uint32_t>(S, num_threads, atoi(argv[3]), atoi(argv[4]));
}