#pragma once

#include <exception>
#include "utils/constant.hpp"
#include <boost/program_options.hpp>
#include <biovoltron/file_io/fasta.hpp>
#include <biovoltron/utility/istring.hpp>
#include <biovoltron/algo/sort/kiss1_sorter.hpp>
#include <biovoltron/algo/sort/kiss2_sorter.hpp>
#include <fstream>
#include <ranges>
#include <variant>

#include <spdlog/spdlog.h>
#include <spdlog/stopwatch.h>

#include "utils/io.hpp"

namespace bpo = boost::program_options;

void suffix_sort_main(
  const bpo::variables_map& generic_vm,
  const bpo::variables_map& command_vm
) {
// TODO: remove me when generic sorting and indexing are supported
  if (generic_vm.count("generic")) {
    throw std::invalid_argument("Generic sorting and indexing are currently not supported.");
  }
  auto fa = std::ifstream{command_vm["fasta"].as<std::string>()};

  auto seq = read_sequence(fa);

  std::ranges::transform(seq, seq.begin(), [](auto& c) { return c % 4; });

  auto k = command_vm["kordered"].as<size_t>();
  auto num_threads = generic_vm["num_threads"].as<size_t>();
  std::variant<biovoltron::KISS1Sorter<uint32_t>, biovoltron::KISS2Sorter<uint32_t>> sorter;
  switch(command_vm["sorting-algorithm"].as<kISS::SortingAlgorithm>()) {
    case kISS::SortingAlgorithm::PARALLEL_SORTING:
      sorter = biovoltron::KISS1Sorter<uint32_t>{};
      break;
    case kISS::SortingAlgorithm::PREFIX_DOUBLING:
      sorter = biovoltron::KISS2Sorter<uint32_t>{};
      break;
    default:
      throw std::invalid_argument("Invalid sorting algorithm");
      break;
  }

  auto S = std::visit([&seq, &k,&num_threads] (auto &&sorter) {
    return sorter.prepare_aligned_ref(seq);
  }, sorter);

  seq.clear();
  seq.shrink_to_fit();

  auto sw = spdlog::stopwatch{};
  auto sa = std::visit([&S, &k,&num_threads] (auto &&sorter) {
    return sorter.get_suffix_array_dna(S, k, num_threads);
  }, sorter);
  SPDLOG_INFO("n = {}, k = {}, suffix sorting elapsed {}", S.size(), k, sw);
}

