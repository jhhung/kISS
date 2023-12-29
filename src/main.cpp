#include <iostream>
#include <sstream>
#include <biovoltron/file_io/fasta.hpp>

// Sample FASTA data
auto fasta = R"(>hello-world
ACTG)";

int main() {
  auto rec = biovoltron::FastaRecord<true>{};
  std::istringstream(fasta) >> rec;

  // Print the name and sequence
  std::cout << rec.name << '\n';
  std::cout << rec.seq << '\n';

  return 0;
}
