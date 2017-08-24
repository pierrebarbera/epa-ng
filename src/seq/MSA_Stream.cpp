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
