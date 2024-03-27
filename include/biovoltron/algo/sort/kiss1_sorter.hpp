#include <biovoltron/algo/sort/core/libsais-sort/sais.hpp>

namespace biovoltron {

template<typename size_type = std::uint32_t>
struct KISS1Sorter {
  static auto
  get_sa(std::string_view ref, size_t sort_len = 256u) {
    auto S = sais::vector<uint8_t>{ref.begin(), ref.end()};
    auto SA = sais::vector<size_type>{};

    // TODO: from parameter
    auto num_threads = std::thread::hardware_concurrency();

    sais::suffix_array(S, SA, sort_len, num_threads);
    return std::vector<size_type>(std::make_move_iterator(SA.begin()), std::make_move_iterator(SA.end()));
  }

  static auto
  get_sa(istring_view ref, size_t sort_len = 256u) {
    auto S = sais::vector<uint8_t>{ref.begin(), ref.end()};
    auto SA = sais::vector<size_type>{};

    // TODO: from parameter
    auto num_threads = std::thread::hardware_concurrency();

    sais::suffix_array(S, SA, sort_len, num_threads);

    return std::vector<size_type>(std::make_move_iterator(SA.begin()), std::make_move_iterator(SA.end()));
  }
};

}
