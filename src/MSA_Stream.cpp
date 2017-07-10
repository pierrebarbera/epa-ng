#include "MSA_Stream.hpp"
#include "file_io.hpp"

void MSA_Stream::read_chunk(const size_t number)
{
  prefetch_chunk_.clear();
  
  if (!fptr_) {
    throw std::runtime_error{"fptr was invalid!"};
  }

  int sites = 0;
  int number_left = std::min(number, max_read_ - num_read_);

  if (number_left == 0) {
    return;
  }

  char * sequence = nullptr;
  char * header = nullptr;
  long sequence_length;
  long header_length;
  long sequence_number;


  while (number_left > 0 and pll_fasta_getnext( fptr_.get(), 
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
    
    prefetch_chunk_.append(header, sequence);
    free(sequence);
    free(header);
    sites = sequence_length;
    number_left--;
  }

  num_read_ += prefetch_chunk_.size();
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

  read_chunk(initial_size);

  // prefetcher_ = std::thread(read_chunk, fptr_.get(), initial_size, std::ref(prefetch_chunk_));
}

size_t MSA_Stream::read_next( MSA_Stream::container_type& result, 
                              const size_t number)
{

#ifdef __PREFETCH
  // join prefetching thread to ensure new chunk exists
  if (prefetcher_.joinable()) {
    prefetcher_.join();
  }
#endif
  // perform pointer swap to data
  std::swap(result, prefetch_chunk_);

  if (!fptr_) {
    throw std::runtime_error{"fptr was invalid during read_next!"};
  }
  
  // start request next chunk from prefetcher (async)
#ifdef __PREFETCH
  prefetcher_ = std::thread(read_chunk, this, number);
#else
  read_chunk(number);
#endif
  // return size of current buffer
  return result.size();
}

MSA_Stream::~MSA_Stream() 
{
#ifdef __PREFETCH
  // avoid dangling threads
  if (prefetcher_.joinable()) {
    prefetcher_.join();
  }
#endif
}
