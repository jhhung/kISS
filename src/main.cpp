#include "utils/options.hpp"
#include <exception>

int main(int argc, char **argv) {
  auto [generic_vm, command_vm] = kISS::argparse(argc, argv);

  throw std::logic_error("Not Implemented");
}
