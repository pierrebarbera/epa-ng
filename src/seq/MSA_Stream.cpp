#include "seq/MSA_Stream.hpp"

#include <chrono>
#include <algorithm>

#include "util/logging.hpp"
#include "net/epa_mpi_util.hpp"

template< typename F >
static void read_chunk_loop(  MSA_Stream::file_type& iter,
                              size_t const number,
                              MSA_Stream::container_type& prefetch_buffer,
                              size_t const max_read,
                              size_t& num_read,
                              F append_func)
{
  prefetch_buffer.clear();

  // auto length = info.sites();
  auto const length = iter->length();

  if( length == 0 ) {
    throw std::runtime_error{"inital length 0?!"};
  }

  size_t number_left = std::min(number, max_read - num_read);

  size_t sites_mismatch = 0;

  while (number_left > 0 and iter)
  {
    auto const sequence_length = iter->length();

    // if (!length) length = sequence_length;

    sites_mismatch += abs( length - sequence_length );

    append_func( prefetch_buffer, iter->label(), iter->sites() );

    number_left--;
    ++iter;
  }

  if ( sites_mismatch ) {
    throw std::runtime_error{"MSA file does not contain equal size sequences"};
  }
  num_read += prefetch_buffer.size();
}

/*
  copy a window from one string to the other. Returns true if this substring was valid
 */
static bool copy_window(std::string const& src,
                        std::string& dest,
                        unsigned int const pos,
                        unsigned int const width,
                        unsigned int const min_nongap_sites=10)
{
  assert( src.size() == dest.size() );

  auto src_start = src.begin() + pos;
  auto src_end = src_start + width;

  auto gap_count = std::count( src_start, src_end, '-' );

  if( width - gap_count < min_nongap_sites ) {
    // substring contained less than the minimum allowed number of nongap chars
    // LOG_DBG << "Skipped pos: " << pos << " width: " << width;
    return false;
  }

  // std:copy( src_start, src_end, dest.begin() + pos );
  dest.replace( pos, width,
                src, pos, width);
  // LOG_DBG << "Multiplexed pos: " << pos << " width: " << width;

  return true;
}

static size_t first_nongap_pos( std::string const& seq )
{
  size_t i = 0;
  while( seq[ i++ ] == '-' ) { }

  return i - 1;
}

/*
  outer function that marries the basic reading loop with the correct "append" function
 */
static void read_chunk( MSA_Stream::file_type& iter,
                        MSA_Info const& info,
                        bool const premasking,
                        size_t const number,
                        MSA_Stream::container_type& prefetch_buffer,
                        size_t const max_read,
                        size_t& num_read,
                        MSA_Stream::multiplex_t const& multiplex)
{
  if( premasking and not multiplex.enable ) {
    read_chunk_loop(iter,
                    number,
                    prefetch_buffer,
                    max_read,
                    num_read,
                    [&gap_mask = info.gap_mask()](MSA_Stream::container_type& buffer,
                                                  std::string const& header,
                                                  std::string const& sequence ) {
                      buffer.append( header, subset_sequence(sequence, gap_mask) );
                    });
  } else if( premasking and multiplex.enable ) {
    read_chunk_loop(iter,
                    number,
                    prefetch_buffer,
                    max_read,
                    num_read,
                    [ &gap_mask = info.gap_mask(),
                      &multiplex](MSA_Stream::container_type& buffer,
                                  std::string const& header,
                                  std::string const& sequence ) {
                      // first get the subsetted sequence on which we operate
                      auto const real_sequence = subset_sequence( sequence, gap_mask );
                      // make a second sequence, equal length, that is all gaps
                      auto multiplex_seq = std::string( real_sequence.size(), '-' );

                      auto const k = multiplex.k;
                      auto const step = multiplex.step;

                      // copy in all possible nonempty kmers
                      for(  size_t i = first_nongap_pos( real_sequence );
                            i < real_sequence.size() - k;
                            i += step ) {
                        // LOG_DBG << "kmer pos: " << i << " width: " << k;
                        
                        if( copy_window( real_sequence, multiplex_seq, i, k ) ) {
                          buffer.append( header, multiplex_seq );
                          
                        }
                        // make sure to reset step worth of chars to gap
                        multiplex_seq.replace( i, step, 1, '-');

                      }
                      
                    });
  } else {
    read_chunk_loop(iter,
                    number,
                    prefetch_buffer,
                    max_read,
                    num_read,
                    []( MSA_Stream::container_type& buffer,
                        std::string const& header,
                        std::string const& sequence ) {
                      buffer.append( header, sequence );
                    });
  }
}

MSA_Stream::MSA_Stream( std::string const& msa_file,
                        MSA_Info const& info,
                        bool const premasking,
                        bool const split)
  : info_(info)
  , premasking_(premasking)
{
  genesis::sequence::FastaReader reader_settings;
  // ensure sequences are uniformly upper case
  reader_settings.site_casing( genesis::sequence::FastaReader::SiteCasing::kToUpper );
  iter_ = genesis::sequence::FastaInputIterator( genesis::utils::from_file( msa_file ), reader_settings );

  if (!iter_) {
    throw std::runtime_error{std::string("Cannot open file: ") + msa_file};
  }

  // if we are under MPI, skip to this ranks assigned part of the input file

  #ifdef __MPI
  if ( split ) {
    // get info about to which sequence to skip to and how much this rank should read
    std::tie(local_seq_offset_, max_read_) = local_seq_package( info.sequences() );

    skip_to_sequence( local_seq_offset_ );
  }
  #else
  static_cast<void>(split);
  #endif

}

size_t MSA_Stream::read_next( MSA_Stream::container_type& result,
                              size_t const number)
{
  if (first_) {//...this is the first chunk
    // then read the first chunk and kick off the next asynchronously
    read_chunk(iter_, info_, premasking_, number, prefetch_chunk_, max_read_, num_read_, multiplex_);
    first_ = false;
  }
#ifdef __PREFETCH
  // join prefetching thread to ensure new chunk exists
  if (prefetcher_.valid()) {
    prefetcher_.wait();
  }
#endif
  // perform pointer swap to data
  std::swap(result, prefetch_chunk_);

  // start request next chunk from prefetcher (async)
#ifdef __PREFETCH
  prefetcher_ = std::async( std::launch::async,
                            read_chunk,
                            std::ref(iter_),
                            std::ref(info_),
                            premasking_,
                            number,
                            std::ref(prefetch_chunk_),
                            max_read_,
                            std::ref(num_read_),
                            std::ref(multiplex_));
#else
  read_chunk(iter_, info_, premasking_, number, prefetch_chunk_, max_read_, num_read_, multiplex_);
#endif
  // return size of current buffer
  return result.size();
}

MSA_Stream::~MSA_Stream()
{
#ifdef __PREFETCH
  // avoid dangling threads
  if (prefetcher_.valid()) {
    prefetcher_.wait();
  }
#endif
}


void MSA_Stream::skip_to_sequence( size_t const n )
{
  // this function is too dirty, disallow usage after first read
  if (not first_) {
    throw std::runtime_error{"Skipping currently not allowed after first read!"};
  }

 #ifdef __PREFETCH
  // join prefetching thread to ensure new chunk exists
  if (prefetcher_.valid()) {
    prefetcher_.wait();
  }
  #endif

  if (n >= num_sequences()) {
    throw std::runtime_error{"Trying to skip out of bounds!"};
  }

  if (n < num_read_) {
    throw std::runtime_error{"Trying to skip behind!"};
  }

  size_t offset = n - num_read_;

  std::advance(iter_, offset);

}
