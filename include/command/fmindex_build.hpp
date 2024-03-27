#pragma once

#include <exception>
#include <boost/program_options.hpp>
#include <biovoltron/algo/align/exact_match/fm_index.hpp>
#include <biovoltron/algo/sort/kpsais_sorter.hpp>

namespace bpo = boost::program_options;

void fmindex_build_main(
  const bpo::variables_map& generic_vm,
  const bpo::variables_map& command_vm
) {
  auto fa_fn = command_vm["fasta"].as<std::string>();
  auto fa = std::ifstream{fa_fn};

  auto seq = biovoltron::istring{};
  auto refs = std::ranges::istream_view<biovoltron::FastaRecord<true>>(fa);
  for (auto &ref : refs)
    seq += ref.seq;

  std::ranges::transform(seq, seq.begin(), [](auto& c) { return c % 4; });

  auto fmi = biovoltron::FMIndex<4, uint32_t, biovoltron::KPsaisSorter<uint32_t>>{
    .LOOKUP_LEN = 0
  };

  fmi.build(seq);

  auto fmi_fout = std::ofstream{fa_fn + ".fmi"};
  fmi.save(fmi_fout);
}

