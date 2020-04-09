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
#include "util/Options.hpp"
#include "io/msa_reader_interface.hpp"

#include "genesis/sequence/formats/fasta_input_iterator.hpp"

class MSA_Stream : public msa_reader
{
public:
  using container_type  = MSA;
  using file_type       = genesis::sequence::FastaInputIterator;

  MSA_Stream (std::string const& msa_file,
              MSA_Info const& info,
              bool const premasking = true,
              bool const split = false);
  MSA_Stream() = default;
  ~MSA_Stream();

  MSA_Stream(MSA_Stream const& other) = delete;
  MSA_Stream(MSA_Stream&& other) = default;

  MSA_Stream& operator= (MSA_Stream const& other) = delete;
  MSA_Stream& operator= (MSA_Stream && other) = default;

  size_t read_next(container_type& result, size_t const number) override;
  size_t num_sequences() const override { return info_.sequences(); }
  size_t local_seq_offset() const override { return local_seq_offset_; }

  void set_multiplex( bool const enable,
                      unsigned int const k = 100,
                      unsigned int const step = 1 )
  {
    multiplex_ = multiplex_t( enable, k, step );
  }

  struct multiplex_t {
    multiplex_t() = default;
    multiplex_t(bool const enable,
                unsigned int const k,
                unsigned int const step)
      : enable( enable )
      , k( k )
      , step( step ) 
      { }

    bool enable = false;
    unsigned int k = 100;
    unsigned int step = 1;
  };

private:
  void skip_to_sequence(size_t const n);

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

  multiplex_t multiplex_;
};
