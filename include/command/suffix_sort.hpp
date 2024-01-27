#pragma once

#include <exception>
#include "utils/constant.hpp"
#include <boost/program_options.hpp>
#include <biovoltron/file_io/fasta.hpp>
#include <biovoltron/utility/istring.hpp>
#include <biovoltron/algo/sort/kpsais_sorter.hpp>
#include <biovoltron/algo/sort/kiss_new_sorter.hpp>
#include <fstream>
#include <ranges>
#include <variant>

#include <spdlog/spdlog.h>
#include <spdlog/stopwatch.h>

namespace bpo = boost::program_options;

void suffix_sort_main(
  const bpo::variables_map& generic_vm,
  const bpo::variables_map& command_vm
) {
  auto fa = std::ifstream{command_vm["fasta"].as<std::string>()};

  auto seq = biovoltron::istring{};
  auto refs = std::ranges::istream_view<biovoltron::FastaRecord<true>>(fa);
  for (auto &ref : refs)
    seq += ref.seq;

  std::ranges::transform(seq, seq.begin(), [](auto& c) { return c % 4; });

  auto k = command_vm["kordered"].as<size_t>();
  auto sw = spdlog::stopwatch{};
  std::variant<biovoltron::KPsaisSorter<uint32_t>, biovoltron::KissNewSorter<uint32_t>> sorter;
  switch(command_vm["sorting-algorithm"].as<kISS::SortingAlgorithm>()) {
    case kISS::SortingAlgorithm::PARALLEL_SORTING:
      sorter = biovoltron::KPsaisSorter<uint32_t>{};
      break;
    case kISS::SortingAlgorithm::PREFIX_DOUBLING:
      sorter = biovoltron::KissNewSorter<uint32_t>{};
      break;
    default:
      throw std::invalid_argument("Invalid sorting algorithm");
      break;
  }
  auto sa = std::visit([&seq, &k] (auto &&sorter) { return sorter.get_sa(seq, k); }, sorter);
  SPDLOG_INFO("k = {}, suffix sorting elapsed {}", k, sw);
}

