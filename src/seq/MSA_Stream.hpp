#pragma once

#include <string>
#include <vector>
#include <stdexcept>
#include <memory>
#include <limits>

#ifdef __PREFETCH
#include <future>
#endif

#include "seq/MSA.hpp"
#include "seq/MSA_Info.hpp"
#include "io/msa_reader_interface.hpp"

#include "genesis/sequence/formats/fasta_input_iterator.hpp"

class MSA_Stream : public msa_reader
{
public:
  using container_type  = MSA;
  using file_type       = genesis::sequence::FastaInputIterator;

  MSA_Stream (const std::string& msa_file,
              const MSA_Info& info,
              const bool premasking = true,
              const bool split = false);
  MSA_Stream() = default;
  ~MSA_Stream();

  MSA_Stream(MSA_Stream const& other) = delete;
  MSA_Stream(MSA_Stream&& other) = default;

  MSA_Stream& operator= (MSA_Stream const& other) = delete;
  MSA_Stream& operator= (MSA_Stream && other) = default;

  size_t read_next(container_type& result, const size_t number) override;
  size_t num_sequences() const override { return info_.sequences(); }
  size_t local_seq_offset() const override { return local_seq_offset_; }

private:
  void skip_to_sequence(const size_t);

private:
  MSA_Info info_;
  file_type iter_;
  // container_type active_chunk_;
  container_type prefetch_chunk_;
#ifdef __PREFETCH
  std::future<void> prefetcher_;
#endif
  bool premasking_ = true;
  size_t num_read_ = 0;
  size_t max_read_ = std::numeric_limits<size_t>::max();
  size_t local_seq_offset_ = 0;
  bool first_ = true;
};
