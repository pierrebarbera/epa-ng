#pragma once

#include <string>
#include <limits>
#include <memory>

#include "MSA.hpp"
#include "encoding.hpp"
#include "template_magic.hpp"

#include "genesis/utils/io/serializer.hpp"
#include "genesis/sequence/formats/fasta_input_iterator.hpp"

constexpr char MAGIC[] = "BFAST\0";
constexpr size_t MAGIC_SIZE = array_size(MAGIC);

using namespace genesis;

class Binary_Fasta
{
private:
  static inline size_t data_section_offset(const size_t num_sequences)
  {
    return MAGIC_SIZE
        + sizeof(num_sequences)
        + (num_sequences *
          sizeof(uint64_t) * 2 );
  }

  static inline void write_header(utils::Serializer& ser, 
                                  const std::vector<size_t>& entry_sizes)
  {
    const uint64_t num_sequences = entry_sizes.size();

    // First part:
    // <magic string><num_sequences>
    ser.put_raw(MAGIC, MAGIC_SIZE);
    ser.put_int(num_sequences);

    // Second part: the offset table for random access
    // <<seqID><seqOffset>><...>
    uint64_t offset = data_section_offset(num_sequences);
    for (uint64_t i = 0; i < num_sequences; ++i) {
      ser.put_int(i);
      ser.put_int(offset);

      // precompute offset of next entry
      const uint64_t sizeof_entry = sizeof(uint64_t)*2 //header_length and seq_length 
                                  + static_cast<uint64_t>(entry_sizes[i]);
      offset += sizeof_entry;
    }
  }

  inline std::vector<size_t> get_entry_sizes(const MSA& msa)
  {
    std::vector<size_t> res;
    for (const auto& s : msa) {
      res.push_back(
        code_.packed_size(s.sequence().size()) +
        s.header().size()
      );
    }
    return res;
  }

  inline void put_encoded(utils::Serializer& ser, 
                          const std::string& seq)
  {
    const auto encoded_seq = code_.to_fourbit(seq);
    // put the size of actual characters
    ser.put_int<uint64_t>(seq.size());
    // pack characters into encoding, write them out
    ser.put_raw_string(encoded_seq);
  }

  std::string get_decoded(utils::Deserializer& des)
  {
    std::string res;

    // get the size of characters that were packed
    const auto decoded_size = des.get_int<uint64_t>();

    // figure out how much that is in bytes
    // (decoded_size = 3 would mean 2 bytes, one for the first two, one for the third plus padding)
    const auto coded_size = code_.packed_size(decoded_size);

    // get the bytes
    auto coded_str = des.get_raw_string(coded_size);

    // decode and return
    return code_.from_fourbit(coded_str, decoded_size);
  }

public:
  Binary_Fasta() = default;
  ~Binary_Fasta() = default;

  void save(const MSA& msa, const std::string& file_name)
  {
    utils::Serializer ser(file_name);

    // Write the header
    auto sizes = get_entry_sizes(msa);
    write_header(ser, sizes);

    // Write the data. Every entry:
    // <header_length (bytes/chars)><header string><sequence_length><encoded sequence padded to next byte>
    // (note: the sequence_length is in number of encoded characters, so for 
    // 4bit the number of bytes to be read is sequence_length * 2) 
    for (const auto& s : msa) {
      ser.put_string(s.header());
      put_encoded(ser, s.sequence());
    }

  }

  MSA load( const std::string& file_name, 
            const size_t skip = 0,
            const size_t number = std::numeric_limits<size_t>::max())
  {
    MSA msa;

    utils::Deserializer des(file_name);

    // read the header info
    char magic[MAGIC_SIZE];
    des.get_raw(magic, MAGIC_SIZE);

    if (strcmp(magic, MAGIC)) {
      throw std::runtime_error{std::string("File is not an epa::Binary_Fasta file: ") + file_name};
    }

    const uint64_t num_sequences = des.get_int<uint64_t>();

    if (skip > num_sequences) {
      throw std::runtime_error{
        std::string("File contains less sequences than you are trying to skip: ") 
        + std::to_string(num_sequences)
        + " vs. "
        + std::to_string(skip)
      };
    }

    std::vector<uint64_t> offset(num_sequences);
    // read the random access table
    for (size_t i = 0; i < num_sequences; ++i) {
      // sequence id
      const auto idx = des.get_int<uint64_t>();
      // offset
      offset[idx] = des.get_int<uint64_t>();
    }

    // skip ahead to correct file offset
    if (skip) {
      // subtract the offset created by the header as skip is
      // essentially a seek from the current position,
      // but the offsets are absolute
      des.skip(offset[skip] - data_section_offset(num_sequences));
    }

    // const size_t to_read = std::min(num_sequences - skip, number);

    for (size_t i = 0; i < number and not des.finished(); ++i) {
      auto label = des.get_string();
      auto sequence = get_decoded(des);
      msa.append( label, sequence );
    }

    return msa;
  }

  void fasta_to_bfast( const std::string& fasta_file, 
                              std::string bfast_file = std::string(""))
  {
    if (!bfast_file.size()) {
      bfast_file = fasta_file + ".bin";
    }

    // detect number of sequences in fasta file
    utils::InputStream instr( std::make_unique< utils::FileInputSource >( fasta_file ));
    auto it = sequence::FastaInputIterator( instr );

    std::vector<size_t> entry_sizes;
    while ( it ) {
      entry_sizes.push_back(
        it->label().size() +
        code_.packed_size(it->sites().size())
        );
      ++it;
    }

    // write the header
    utils::Serializer ser(bfast_file);
    write_header(ser, entry_sizes);

    // write the data
    utils::InputStream instr_again( std::make_unique< utils::FileInputSource >( fasta_file ));
    auto it_again = sequence::FastaInputIterator( instr_again );
    while ( it_again ) {
      ser.put_string(it_again->label());
      put_encoded(ser, it_again->sites());
      ++it_again;
    }
  }

private:
  FourBit code_;
  
};
