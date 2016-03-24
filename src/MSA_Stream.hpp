#pragma once

#include <string>
#include <vector>

#include "pllhead.hpp"
#include "Sequence.hpp"

class MSA_Stream {
public:
  MSA_Stream (std::string& msa_file)
  {
    if(!(fptr_ = pll_fasta_open(msa_file.c_str(), pll_map_fasta)))
      throw runtime_error{"Couldn't open MSA file."};
  }
  ~MSA_Stream ()
  {
    pll_fasta_close(fptr_);
  }
  void clear() {sequence_list_.clear();}

  unsigned int read_next(const unsigned int number)
  {
    int sites = 0;
    unsigned int number_left = number;

    while (number_left and pll_fasta_getnext(file, &header, &header_length, &sequence, &sequence_length, &sequence_number))
    {
      if (sites and sites not sequence_length)
        throw runtime_error{"MSA file does not contain equal size sequences"};

      if (!sites) sites = sequence_length;

      sequence_list_.emplace_back(header, sequence);
      free(sequence);
      free(header);
      sites = sequence_length;
      number_left--;
    }
    return number - number_left;
  }

  const Sequence& operator[](const unsigned int idx) const
  {
    if((size_t) idx >= sequence_list_.size())
      throw runtime_error{string("Trying to access MSA entry out of bounds. i = ") + to_string(i) };

    return sequence_list_[idx];
  }

private:
  FILE* fptr_;
  std::vector<Sequence> sequence_list_;

};
