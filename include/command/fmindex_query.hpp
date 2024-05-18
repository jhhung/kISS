#pragma once

#include <sstream>
#include <exception>
#include <boost/program_options.hpp>

#include "utils/io.hpp"

namespace bpo = boost::program_options;

void fmindex_query_main(
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

  auto fmi_fin = std::ifstream{fa_fn + ".fmi"};
  assert(fmi_fin);
  fmi.load(fmi_fin);

  if (command_vm.count("query")) {
    auto query = biovoltron::Codec::to_istring(command_vm["query"].as<std::string>());

    auto [beg, end, offs] = fmi.get_range(query);

    auto positions = fmi.get_offsets(beg, end);
    auto headn = command_vm["headn"].as<size_t>();

    auto ending = [](int x) {
      x %= 100;
      if (x / 10 == 1)
        return "th";
      if (x % 10 == 1)
        return "st";
      if (x % 10 == 2)
        return "nd";
      if (x % 10 == 3)
        return "rd";
      return "th";
    };

    SPDLOG_INFO("query = {} found {} times", biovoltron::Codec::to_string(query), positions.size());
    for (auto i = size_t{}; i < std::min(headn, positions.size()); i++) {
      auto loc = positions[i];
      SPDLOG_INFO(
        "The {}-{} position is {}, content of substring is {}",
        i + 1, ending(i + 1), loc,
        biovoltron::Codec::to_string(seq.substr(loc, query.size()))
      );
    }
  }

  if (command_vm.count("batch")) {
    auto pfile = std::ifstream{command_vm["batch"].as<std::string>()};
    uint32_t query_len;
    uint32_t num_query;
    pfile.read(reinterpret_cast<char*>(&query_len), sizeof(uint32_t));
    pfile.read(reinterpret_cast<char*>(&num_query), sizeof(uint32_t));

    SPDLOG_INFO("query_len: {}, num_query: {}", query_len, num_query);
    auto query = std::string{};
    auto time = double{};
    auto occ = size_t{};
    auto location_checksum = size_t{};
    query.resize(query_len);
    while (num_query--) {
      pfile.read(query.data(), query_len);
      auto iq = biovoltron::Codec::to_istring(query);
      auto tic = std::chrono::high_resolution_clock::now();

      auto [beg, end, offs] = fmi.get_range(iq);
      auto positions = fmi.get_offsets(beg, end);
      auto toc = std::chrono::high_resolution_clock::now();
      occ += positions.size();
      time += (toc - tic).count() / 1e9;

      if (num_query % 100000 == 0)
        SPDLOG_DEBUG("remain: {}, time: {}", num_query, time);

      for (const auto &v : positions)
        location_checksum += v;
    }
    SPDLOG_INFO("searching time: {} seconds", time);
    SPDLOG_INFO("number of matched locations: {}", occ);
    SPDLOG_INFO("location checksum: {}", location_checksum);
  }
}

