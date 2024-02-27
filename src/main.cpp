// SPDLOG LOGGER
#define SPDLOG_ACTIVE_LEVEL SPDLOG_LEVEL_DEBUG

#include "command/suffix_sort.hpp"
#include "command/fmindex_build.hpp"
#include "command/fmindex_query.hpp"

#include "utils/options.hpp"
#include <exception>

#include <iostream>
#include <string>
#include <map>

namespace bpo = boost::program_options;

int main(int argc, char **argv) {
  auto [generic_vm, command_vm] = kISS::argparse(argc, argv);

  auto command = generic_vm["command"].as<std::string>();

  auto dispatcher = std::map<
    std::string,
    std::function<void(bpo::variables_map&, bpo::variables_map&)>
  >{
    {"suffix_sort",   suffix_sort_main},
    {"fmindex_build", fmindex_build_main},
    {"fmindex_query", fmindex_query_main}
  };

  dispatcher[command](generic_vm, command_vm);
}
