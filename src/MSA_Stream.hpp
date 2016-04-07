#pragma once

#include <string>
#include <vector>
#include <stdexcept>
#include <memory>

#include "pllhead.hpp"
#include "Sequence.hpp"
#include "pll_util.hpp"

class MSA_Stream {
public:
  MSA_Stream (const std::string& msa_file);
  MSA_Stream() : fptr_(nullptr, fasta_close) { }
  ~MSA_Stream () = default;

  MSA_Stream(MSA_Stream const& other) = delete;
  MSA_Stream(MSA_Stream&& other) = default;

  MSA_Stream& operator= (MSA_Stream const& other) = delete;
  MSA_Stream& operator= (MSA_Stream && other) = default;

  void clear() { sequence_list_.clear(); }
  unsigned int read_next(const unsigned int number);
  const Sequence& operator[](const unsigned int idx) const;

private:
  std::unique_ptr<pll_fasta_t, fasta_deleter> fptr_;
  std::vector<Sequence> sequence_list_;

};
