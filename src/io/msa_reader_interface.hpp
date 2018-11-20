#pragma once

#include "seq/MSA.hpp"

class msa_reader
{

public:
  msa_reader() = default;
  virtual ~msa_reader() = default;

  virtual size_t num_sequences() const = 0;
  virtual size_t local_seq_offset() const = 0;
  virtual size_t read_next(MSA& result, const size_t number) = 0;

};
