#pragma once

// for TBB deprecated warning
#define TBB_SUPPRESS_DEPRECATED_MESSAGES 1

// SPDLOG LOGGER
#define SPDLOG_ACTIVE_LEVEL SPDLOG_LEVEL_DEBUG

#include "constant.hpp"
#include <boost/program_options.hpp>

#include "spdlog/spdlog.h"
#include "spdlog/sinks/stdout_color_sinks.h"

#include <iostream>
#include <thread>

namespace kISS {

namespace bpo = boost::program_options;

auto generic_options = []{ // {{{
  auto generic_options = bpo::options_description{"Generic options"};
  generic_options.add_options()(
    "help,h",
    "produce help message"
  )(
    "version,v",
    "print version string"
  )(
    "num_threads,t",
    bpo::value<size_t>()
      ->default_value(std::thread::hardware_concurrency())
      ->value_name("NUM"),
    "number of thread"
  )(
    "verbose",
    "print more information"
  );
  return generic_options;
}(); // }}}
auto [generic_options_cmdline, generic_positional] = []{ // {{{
  auto generic_options_hidden = bpo::options_description{};
  generic_options_hidden.add_options()(
    "command",
    bpo::value<std::string>(),
    "Command to execute"
  )(
    "subargs",
    bpo::value<std::vector<std::string>>(),
    "Arguments for command"
  );

  auto generic_positional = bpo::positional_options_description{};
  generic_positional.add("command", 1).add("subargs", -1);

  auto generic_options_cmdline = bpo::options_description{};
  generic_options_cmdline.add(generic_options).add(generic_options_hidden);

  return std::tuple{generic_options_cmdline, generic_positional};
}(); // }}}

auto suffix_sort_options = []{ // {{{
  namespace bpo = boost::program_options;
  auto suffix_sort_options = bpo::options_description{
    "\n./kISS suffix_sort [--option ...] <FASTA/text filename>\n\n"
    "Options"
  };
  suffix_sort_options.add_options()
  (
    ",k",
    bpo::value<size_t>()
      ->value_name("NUM")
      ->default_value(256),
    "a k-ordered value, where each suffix is sorted based on the first k characters.\n"
    "Using -1 indicates unbounded sorting."
  );
  return suffix_sort_options;
}(); // }}}
auto [suffix_sort_options_cmdline, suffix_sort_positional] = []{ // {{{
  auto suffix_sort_options_hidden = bpo::options_description{};
  suffix_sort_options_hidden.add_options()(
    "fasta",
    bpo::value<std::string>()->required(),
    "fasta file path"
  );

  auto suffix_sort_positional = bpo::positional_options_description{};
  suffix_sort_positional.add("fasta", 1);

  auto suffix_sort_options_cmdline = bpo::options_description{};
  suffix_sort_options_cmdline.add(suffix_sort_options).add(suffix_sort_options_hidden);

  return std::tuple{suffix_sort_options_cmdline, suffix_sort_positional};
}(); // }}}

auto fmindex_build_options = []{ // {{{
  namespace bpo = boost::program_options;
  auto fmindex_build_options = bpo::options_description{
    "\n./kISS fmindex_build [--option ...] <FASTA/text filename>\n\n"
    "Options"
  };
  fmindex_build_options.add_options()
  (
    ",k",
    bpo::value<size_t>()
      ->value_name("NUM")
      ->default_value(256),
    "Sets the maximum query length for the fmindex search as k - sa_sample_rate; using -1 indicates an unlimited query length."
  );
  return fmindex_build_options;
}(); // }}}
auto [fmindex_build_options_cmdline, fmindex_build_positional] = []{ // {{{
  auto fmindex_build_options_hidden = bpo::options_description{};
  fmindex_build_options_hidden.add_options()(
    "fasta",
    bpo::value<std::string>()->required(),
    "fasta file path"
  );

  auto fmindex_build_positional = bpo::positional_options_description{};
  fmindex_build_positional.add("fasta", 1);

  auto fmindex_build_options_cmdline = bpo::options_description{};
  fmindex_build_options_cmdline.add(fmindex_build_options).add(fmindex_build_options_hidden);

  return std::tuple{fmindex_build_options_cmdline, fmindex_build_positional};
}(); // }}}

auto fmindex_query_options = []{ // {{{
  namespace bpo = boost::program_options;
  auto fmindex_query_options = bpo::options_description{
    "./kISS fmindex_query [--option ...] <FASTA/text filename>\n\n"
    "Options"
  };
  fmindex_query_options.add_options()
  (
    ",q",
    bpo::value<std::string>()
      ->value_name("STR"),
    "Content of the query string."
  )
  (
    ",f",
    bpo::value<std::string>()
      ->value_name("STR"),
    "File containing multiple query string contents, with each query string separated by a newline."
  );
  return fmindex_query_options;
}(); // }}}
auto [fmindex_query_options_cmdline, fmindex_query_positional] = []{ // {{{
  auto fmindex_query_options_hidden = bpo::options_description{};
  fmindex_query_options_hidden.add_options()(
    "fasta",
    bpo::value<std::string>()->required(),
    "fasta file path"
  );

  auto fmindex_query_positional = bpo::positional_options_description{};
  fmindex_query_positional.add("fasta", 1);

  auto fmindex_query_options_cmdline = bpo::options_description{};
  fmindex_query_options_cmdline.add(fmindex_query_options).add(fmindex_query_options_hidden);

  return std::tuple{fmindex_query_options_cmdline, fmindex_query_positional};
}(); // }}}

auto help_options = [] { // {{{
  namespace bpo = boost::program_options;
  auto help_options = bpo::options_description{
    BANNER + '\n' +
    "kISS [--generic-option ...] cmd [--cmd-specific-option ...]"
  };
  help_options
    .add(generic_options)
    .add(suffix_sort_options)
    .add(fmindex_build_options)
    .add(fmindex_query_options);
  return help_options;
}(); // }}}

template <typename... Args> // {{{
auto parse(
  bpo::options_description &options,
  bpo::positional_options_description &positional,
  bool allow_unregistered,
  const Args&... args
) {
  try {
    auto parsed = bpo::command_line_parser(args...)
      .options(options)
      .positional(positional);

    if (allow_unregistered)
      parsed.allow_unregistered();

    auto parsed_run = parsed.run();
    auto vm = bpo::variables_map{};
    bpo::store(parsed_run, vm);
    bpo::notify(vm);

    return std::tuple{parsed_run, vm};
  } catch (bpo::error &e) {
    std::cerr << help_options << std::endl;
    std::cerr << e.what() << std::endl;
    exit(1);
  }
} // }}}
auto argparse(int argc, char **argv) { // {{{
  namespace bpo = boost::program_options;

  auto [generic_parsed, generic_vm] = parse(
    generic_options_cmdline,
    generic_positional,
    true, argc, argv
  );

  if (generic_vm.count("version")) {
    std::cerr << VERSION << std::endl;
    exit(1);
  }

  // invalid argument or help
  if (not generic_vm.count("command") or generic_vm.count("help")) {
    std::cerr << help_options << std::endl;
    exit(1);
  }

  // setup spdlog into stderr and verbosity
  auto err_logger = spdlog::stderr_color_mt("stderr");
  spdlog::set_default_logger(err_logger);
  if (generic_vm.count("verbose"))
    spdlog::set_level(spdlog::level::debug);

  auto opts = bpo::collect_unrecognized(
    generic_parsed.options, bpo::include_positional
  );

  // remove "command"
  opts.erase(opts.begin());

  auto command = generic_vm["command"].as<std::string>();
  if (command == "suffix_sort") {
    auto [suffix_sort_parsed, suffix_sort_vm] = parse(
      suffix_sort_options_cmdline, suffix_sort_positional, false, opts
    );
    return std::tuple{generic_vm, suffix_sort_vm};
  } else if (command == "fmindex_build") {
    auto [fmindex_build_parsed, fmindex_build_vm] = parse(
      fmindex_build_options_cmdline, fmindex_build_positional, false, opts
    );
    return std::tuple{generic_vm, fmindex_build_vm};
  } else if (command == "fmindex_query") {
    auto [fmindex_query_parsed, fmindex_query_vm] = parse(
      fmindex_query_options_cmdline, fmindex_query_positional, false, opts
    );
    return std::tuple{generic_vm, fmindex_query_vm};
  }

  // TODO: help message per command
  std::cerr << help_options << std::endl;
  throw bpo::invalid_option_value(command);
} // }}}

}
