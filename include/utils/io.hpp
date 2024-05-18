#pragma once

#include <string>
#include <biovoltron/utility/istring.hpp>

auto read_sequence(auto &fin) {
  auto seq = biovoltron::istring{};
  if (fin.peek() == '>') { // fasta mode
    auto refs = std::ranges::istream_view<biovoltron::FastaRecord<true>>(fin);
    for (auto &ref : refs)
      seq += ref.seq;
  } else { // text mode
    auto line = std::string{};
    while (getline(fin, line))
      seq += biovoltron::Codec::to_istring(line);
  }
  return seq;
}
