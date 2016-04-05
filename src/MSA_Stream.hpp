#pragma once

#include <string>
#include <vector>
#include <stdexcept>

#include "pllhead.hpp"
#include "Sequence.hpp"

class MSA_Stream {
public:
  MSA_Stream (const std::string& msa_file);
  MSA_Stream() : fptr_(nullptr) { }
  ~MSA_Stream () { pll_fasta_close(fptr_); }

  void clear() { sequence_list_.clear(); }
  unsigned int read_next(const unsigned int number);
  const Sequence& operator[](const unsigned int idx) const;

private:
  pll_fasta_t * fptr_;
  std::vector<Sequence> sequence_list_;

};
