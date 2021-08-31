#pragma once

#include <biovoltron/utility/istring.hpp>

namespace biovoltron {

template<bool Encoded = false>
struct FastaRecord {
  constexpr static auto encoded = Encoded;
  constexpr static auto START_SYMBOL = '>';
  std::string name;
  std::conditional_t<Encoded, istring, std::string> seq;

  operator auto() const {
    if constexpr (Encoded)
      return FastaRecord<!Encoded>{name, Codec::to_string(seq)};
    else
      return FastaRecord<!Encoded>{name, Codec::to_istring(seq)};
  }
};

template<bool Encoded = false>
struct FastqRecord : FastaRecord<Encoded> {
  constexpr static auto START_SYMBOL = '@';
  constexpr static auto DELIM = '+';
  std::string qual;
};

template<class R>
  requires std::derived_from<R, FastaRecord<R::encoded>>
inline auto&
operator>>(std::istream& is, R& record) {
  if (is >> std::ws; is.peek() != record.START_SYMBOL) {
    is.clear(std::ios::failbit);
    return is;
  }

  auto line = std::string{};
  std::getline(is, line);
  record.name = line.substr(1, line.find_first_of(" \t", 1) - 1);
  for (record.seq.clear(); std::getline(is, line);) {
    if constexpr (R::encoded)
      record.seq += Codec::to_istring(line);
    else
      record.seq += line;
    if constexpr (std::same_as<R, FastqRecord<R::encoded>>) {
      if (is.peek() == record.DELIM)
        break;
    } else {
      if (is.peek() == record.START_SYMBOL)
        return is;
    }
  }
  if constexpr (std::same_as<R, FastqRecord<R::encoded>>) {
    getline(is, line);
    for (record.qual.clear(); getline(is, line);) {
      record.qual += line;
      if (is.peek() == record.START_SYMBOL)
        return is;
    }
  }
  is.clear();
  return is;
}

template<class R>
  requires std::derived_from<R, FastaRecord<R::encoded>>
inline auto&
operator<<(std::ostream& os, const R& record) {
  os << record.START_SYMBOL << record.name << "\n";
  if constexpr (R::encoded)
    os << Codec::to_string(record.seq);
  else
    os << record.seq;
  if constexpr (std::same_as<R, FastqRecord<R::encoded>>)
    os << "\n" << record.DELIM << "\n" << record.qual;
  return os;
}

}  // namespace biovoltron
