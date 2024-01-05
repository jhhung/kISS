#pragma once

#include <exception>
#include <boost/program_options.hpp>

namespace bpo = boost::program_options;

void fmindex_query_main(
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

  auto fmi_fin = std::ifstream{fa_fn + ".fmi"};
  assert(fmi_fin);
  fmi.load(fmi_fin);

  auto query = biovoltron::Codec::to_istring(command_vm["query"].as<std::string>());

  auto [beg, end, offs] = fmi.get_range(query);

  auto positions = fmi.get_offsets(beg, end);

  SPDLOG_INFO("query = {} found {} times", biovoltron::Codec::to_string(query), positions.size());
}

