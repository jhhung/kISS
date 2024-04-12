#pragma once

#include <string>
#include <fstream>
#include <vector>

auto read_fasta(std::string filename) {
  auto fin = std::ifstream(filename);
  auto line = std::string{};
  auto seq = std::string{};
  while (std::getline(fin, line)) {
    if (line[0] == '>')
      continue;
    for (auto &c : line) {
      seq.push_back(std::toupper(c));
    }
  }
  return seq;
}

auto read_fasta_dna(std::string filename) {
  auto fin = std::ifstream(filename);
  auto line = std::string{};
  auto seq = std::string{};
  while (std::getline(fin, line)) {
    if (line[0] == '>')
      continue;
    for (auto &c : line) {
      c = std::toupper(c);
      if (c == 'A')
        c = 1;
      else if (c == 'C')
        c = 2;
      else if (c == 'G')
        c = 3;
      else if (c == 'T')
        c = 4;
      else
        c = 1;
      seq.push_back(std::toupper(c));
    }
  }
  return seq;
}