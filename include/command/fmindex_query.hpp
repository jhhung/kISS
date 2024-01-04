#pragma once

#include <exception>
#include <boost/program_options.hpp>

namespace bpo = boost::program_options;

void fmindex_query_main(
  const bpo::variables_map& generic_vm,
  const bpo::variables_map& command_vm
) {
  throw std::logic_error("Not Implemented");
}

