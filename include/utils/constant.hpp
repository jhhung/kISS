#pragma once

#include <string>
#include <stdexcept>

namespace kISS {

const auto VERSION = std::string{"1.0.0"};

const auto BANNER = std::string{R"( _     ___  ____  ____
| | __|_ _|/ ___|/ ___|
| |/ / | | \___ \\___ \
|   <  | |  ___) |___) |
|_|\_\|___||____/|____/ )" + VERSION + '\n'};

enum SortingAlgorithm
{
  PARALLEL_SORTING,
  PREFIX_DOUBLING
};

}
