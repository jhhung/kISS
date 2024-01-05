#pragma once

#include <biovoltron/utility/istring.hpp>


namespace biovoltron {

/**
 * @ingroup file_io
 * @brief This is a data structure for fasta file.
 * @tparam Encoded seq type == istring if Encoded else string
 */
template<bool Encoded = false>
struct FastaRecord {
  /**
   * @brief
   * boolean to define encode type
   *
   * - encoded = false: encoded by string (ACGT)
   * - encoded = true: encoded by integer (0123)
   */
  constexpr static auto encoded = Encoded;

  /**
   * @brief
   * Start symbol of reference genome
   *
   * - using ">" in fasta file type
   */
  constexpr static auto START_SYMBOL = '>';

  /**
   * @brief
   * Name of reference genome (or reads in FastqRecord)
   */
  std::string name;

  /**
   * @brief
   * Sequence of reference genome (or reads in FastqRecord)
   */
  std::conditional_t<Encoded, istring, std::string> seq;

  /**
   * @brief
   * implicit conversion operator (converted to FastaRecord<!Encoded>)
   *
   * - If you have existing FastaRecord (or FastqRecord) object a = FastaRecord<false>,
   *   and you want to change a's encoding from std::string to biovoltron::istring,
   *   this operator allows you to implicitly convert you object to the !Encoded version.
   * Example
   * ```cpp
   *  FastaRecord<false> a;     // a is encoded by std::string
   *  FastaRecord<true> b = a;  // b is encoded by biovoltron::istring (using implicit conversion)
   * ```
   */
  operator auto() const {
    if constexpr (Encoded)
      return FastaRecord<!Encoded>{name, Codec::to_string(seq)};
    else
      return FastaRecord<!Encoded>{name, Codec::to_istring(seq)};
  }
};

/**
 * @ingroup file_io
 * @brief This is a data structure for fastq file.
 *
 * This is a data structure for fastq file.
 * - inherit from FastaRecord
 * - for fastq file type to record reads after sequencing.
 *   - include reads and reads quality.
 */
template<bool Encoded = false>
struct FastqRecord : FastaRecord<Encoded> {
  /**
   * @brief
   * Start symbol of reads
   *
   * - using "@" in fastq file type
   */
  constexpr static auto START_SYMBOL = '@';

  /**
   * @brief
   * Delimiter between reads and quality of reads
   *
   * - using "+" in fastq file type
   */
  constexpr static auto DELIM = '+';

  /**
   * @brief
   * Quality of reads
   *
   * - Quality record by character from "!" to "~"
   *   - encoded by ASCII from 33 to 126 transfer to quality 0~93
   *     - example: "!" ASCIII = (0+33), quality = 0; "?" ASCIII = (30+33), quality = 30
   *   - Q = -10*log(P)
   *     - Q: quality
   *     - P: error probability
   */
  std::string qual;
};


/**
 * @brief
 * read fasta or fastq file and record data into FastaRecord or FastqRecord
 *
 * - class R "record": must be FastaRecord or FastqRecord (inherit from FastaRecord)
 * - every single genome(reads) will record one FastaRecord(FastqRecord)
 * - All sequences will record into record.seq
 * - if FastqRecord type, program will check record.DELIM and record quality of reads into record.qual
 */
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


/**
 * @brief
 * output sequence and quality data in FastaRecord or FastqRecord
 *
 * - class R "record": must be FastaRecord or FastqRecord (inherit from FastaRecord)
 * - All sequences in record.seq will output
 * - if FastqRecord type, program will check record.DELIM and output quality of reads in record.qual.
 */
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
