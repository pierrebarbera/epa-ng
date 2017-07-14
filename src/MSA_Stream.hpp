#pragma once

#include <string>
#include <vector>
#include <stdexcept>
#include <memory>
#include <limits>

#ifdef __PREFETCH
#include <future>
#endif

#include "pllhead.hpp"
#include "Sequence.hpp"
#include "pll_util.hpp"
#include "MSA.hpp"

class MSA_Stream {
public:
  using container_type  = MSA;
  using file_type       = std::unique_ptr<pll_fasta_t, fasta_deleter>;

  MSA_Stream (const std::string& msa_file, 
              const size_t initial_size,
              const size_t offset=0,
              const size_t max_read=std::numeric_limits<size_t>::max());
  MSA_Stream() : fptr_(nullptr, fasta_close) { }
  ~MSA_Stream();

  MSA_Stream(MSA_Stream const& other) = delete;
  MSA_Stream(MSA_Stream&& other) = default;

  MSA_Stream& operator= (MSA_Stream const& other) = delete;
  MSA_Stream& operator= (MSA_Stream && other) = default;

  size_t read_next(container_type& result, const size_t number);

private:
  file_type fptr_;
  // container_type active_chunk_;
  container_type prefetch_chunk_;
#ifdef __PREFETCH
  std::future<void> prefetcher_;
#endif
  size_t num_read_ = 0;
  size_t max_read_;
};
