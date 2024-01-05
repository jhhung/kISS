#pragma once

#include <exception>
#include <boost/program_options.hpp>
#include <biovoltron/file_io/fasta.hpp>
#include <biovoltron/utility/istring.hpp>
#include <biovoltron/algo/sort/kpsais_sorter.hpp>
#include <fstream>
#include <ranges>

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

  auto sorter = biovoltron::KPsaisSorter<uint32_t>{};

  auto k = command_vm["kordered"].as<size_t>();
  auto sw = spdlog::stopwatch{};
  auto sa = sorter.get_sa(seq, k);
  SPDLOG_INFO("k = {}, suffix sorting elapsed {}", k, sw);
}

