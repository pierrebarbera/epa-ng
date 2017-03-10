#pragma once

#include <string>
#include <vector>
#include <stdexcept>
#include <memory>

#ifdef __PREFETCH
#include <thread>
#endif

#include "pllhead.hpp"
#include "Sequence.hpp"
#include "pll_util.hpp"

class MSA_Stream {
public:
  using container_type  = std::vector<Sequence>;
  using file_type       = std::unique_ptr<pll_fasta_t, fasta_deleter>;

  MSA_Stream (const std::string& msa_file , const size_t initial_size);
  MSA_Stream() : fptr_(nullptr, fasta_close) { }
  ~MSA_Stream();

  MSA_Stream(MSA_Stream const& other) = delete;
  MSA_Stream(MSA_Stream&& other) = default;

  MSA_Stream& operator= (MSA_Stream const& other) = delete;
  MSA_Stream& operator= (MSA_Stream && other) = default;

  // void clear() { active_chunk_->clear(); }
  size_t read_next(const size_t number);
  const Sequence& operator[](const size_t idx) const;

private:
  file_type fptr_;
  container_type active_chunk_;
  container_type prefetch_chunk_;
#ifdef __PREFETCH
  std::thread prefetcher_;
#endif
};
