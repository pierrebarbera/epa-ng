#include "seq/MSA_Stream.hpp"
#include "io/file_io.hpp"

#include <chrono>

static void read_chunk( MSA_Stream::file_type::pointer fptr, 
                        const size_t number, 
                        MSA_Stream::container_type& prefetch_buffer,
                        const size_t max_read,
                        size_t& num_read)
{
  prefetch_buffer.clear();

  if (!fptr) {
    throw std::runtime_error{"fptr was invalid!"};
  }

  int sites = 0;
  int number_left = std::min(number, max_read - num_read);

  if (number_left <= 0) {
    return;
  }

  char * sequence = nullptr;
  char * header = nullptr;
  long sequence_length;
  long header_length;
  long sequence_number;

  while (number_left > 0 and pll_fasta_getnext( fptr, 
                                            &header, 
                                            &header_length, 
                                            &sequence, 
                                            &sequence_length, 
                                            &sequence_number))
  {
    if (sites and (sites != sequence_length)) {
      throw std::runtime_error{"MSA file does not contain equal size sequences"};
    }

    if (!sites) sites = sequence_length;

    for (long i = 0; i < sequence_length; ++i) {
      sequence[i] = toupper(sequence[i]);
    }
    
    prefetch_buffer.append(header, sequence);
    free(sequence);
    free(header);
    sites = sequence_length;
    number_left--;
  }

  num_read += prefetch_buffer.size();
}

MSA_Stream::MSA_Stream( const std::string& msa_file, 
                        const size_t initial_size,
                        const size_t offset,
                        const size_t max_read)
  : fptr_(nullptr, fasta_close)
  , max_read_(max_read)
  , filename_(msa_file)
  , initial_size_(initial_size)
{
  fptr_ = file_type(pll_fasta_open(msa_file.c_str(), pll_map_fasta),
                    fasta_close);
  if (!fptr_) {
    throw std::runtime_error{std::string("Cannot open file: ") + msa_file};
  }

  if (offset) {
    if (pll_fasta_fseek(fptr_.get(), offset, SEEK_SET)) {
      throw std::runtime_error{"Unable to fseek on the fasta file."};
    }
  }

  read_chunk(fptr_.get(), initial_size, prefetch_chunk_, max_read_, num_read_);

  // prefetcher_ = std::thread(read_chunk, fptr_.get(), initial_size, std::ref(prefetch_chunk_));
}

size_t MSA_Stream::read_next( MSA_Stream::container_type& result, 
                              const size_t number)
{

#ifdef __PREFETCH
  // join prefetching thread to ensure new chunk exists
  if (prefetcher_.valid()) {
    prefetcher_.wait();
  }
#endif
  // perform pointer swap to data
  std::swap(result, prefetch_chunk_);

  if (!fptr_) {
    throw std::runtime_error{"fptr was invalid during read_next!"};
  }
  
  // start request next chunk from prefetcher (async)
#ifdef __PREFETCH
  prefetcher_ = std::async( std::launch::async,
                            read_chunk, 
                            fptr_.get(), 
                            number, 
                            std::ref(prefetch_chunk_),
                            max_read_,
                            std::ref(num_read_));
#else
  read_chunk(fptr_.get(), number, prefetch_chunk_, max_read_, num_read_);
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

void MSA_Stream::constrain(const size_t max_read)
{
  max_read_ = max_read;
}

void MSA_Stream::skip_to_sequence(const size_t n)
{
 #ifdef __PREFETCH
  // join prefetching thread to ensure new chunk exists
  if (prefetcher_.valid()) {
    prefetcher_.wait();
  }
#endif 
  // ensure we have the offsets
  if ( offsets_.size()  == 0
    or num_sequences_   == 0 ) {
    num_sequences_ = num_sequences();
  }

  if (n >= offsets_.size()) {
    throw std::runtime_error{"Trying to skip out of bounds!"};
  }
  auto offset = offsets_[n];

  // seek the fileptr
  if (offset) {
    if (pll_fasta_fseek(fptr_.get(), offset, SEEK_SET)) {
      throw std::runtime_error{"Unable to fseek on the fasta file."};
    }
  }

  if (num_read_ == initial_size_) {
    num_read_ = 0;
  }

  // kick off reading a chunk
#ifdef __PREFETCH
  prefetcher_ = std::async( std::launch::async,
                            read_chunk, 
                            fptr_.get(), 
                            initial_size_, 
                            std::ref(prefetch_chunk_),
                            max_read_,
                            std::ref(num_read_));
#else
  read_chunk(fptr_.get(), number, prefetch_chunk_, max_read_, num_read_);
#endif
}

#include <fstream>

size_t MSA_Stream::num_sequences()
{
  assert(filename_.size());
  size_t ret = num_sequences_;

  if (ret == 0) {
    // go through entire file, count sequences ('>' character)
    std::ifstream stream(filename_);
    
    for (size_t pos = 0; stream.peek() != EOF; ++pos) {
      char c;
      stream.get(c);
      // TODO not good enough, need lots of escaping. use existing fasta parser, calculate the offset
      if (c == '>') {
        offsets_.push_back(pos);
      }
    }
    ret = offsets_.size();
  }

  return ret;
}
