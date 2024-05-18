#pragma once

#include <exception>
#include <boost/program_options.hpp>
#include <biovoltron/algo/align/exact_match/fm_index.hpp>
#include <biovoltron/algo/sort/kiss1_sorter.hpp>

#include "utils/io.hpp"

namespace bpo = boost::program_options;

void fmindex_build_main(
  const bpo::variables_map& generic_vm,
  const bpo::variables_map& command_vm
) {
// TODO: remove me when generic sorting and indexing are supported
  if (generic_vm.count("generic")) {
    throw std::invalid_argument("Generic sorting and indexing are currently not supported.");
  }
  auto fa_fn = command_vm["fasta"].as<std::string>();
  auto fa = std::ifstream{fa_fn};

  auto seq = read_sequence(fa);

  std::ranges::transform(seq, seq.begin(), [](auto& c) { return c % 4; });

  auto fmi = biovoltron::FMIndex<4, uint32_t, biovoltron::KISS1Sorter<uint32_t>>{
    .LOOKUP_LEN = 0
  };

  fmi.build(seq);

  auto fmi_fout = std::ofstream{fa_fn + ".fmi"};
  fmi.save(fmi_fout);
}

