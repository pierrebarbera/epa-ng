#include "MSA_Stream.hpp"

using namespace std;

static void read_chunk(MSA_Stream::file_type& fptr, const size_t number, MSA_Stream::container_type& prefetch_buffer)
{
  assert(fptr);
  int sites = 0;
  size_t number_left = number;

  char * sequence = nullptr;
  char * header = nullptr;
  long sequence_length;
  long header_length;
  long sequence_number;

  while (number_left and
    pll_fasta_getnext(fptr.get(), &header, &header_length, &sequence, &sequence_length, &sequence_number))
  {
    if (sites && (sites != sequence_length))
      throw runtime_error{"MSA file does not contain equal size sequences"};

    if (!sites) sites = sequence_length;

    for (long i = 0; i < sequence_length; ++i) sequence[i] = toupper(sequence[i]);
    
    prefetch_buffer.emplace_back(header, sequence);
    free(sequence);
    free(header);
    sites = sequence_length;
    number_left--;
  }
  //return number - number_left; 
}

MSA_Stream::MSA_Stream (const string& msa_file, const size_t initial_size)
  : fptr_(nullptr, fasta_close)
{
  fptr_ = unique_ptr<pll_fasta_t, fasta_deleter>(
                      pll_fasta_open(msa_file.c_str(), pll_map_fasta),
                      fasta_close);
  if(!fptr_)
    throw runtime_error{string("Cannot open file: ") + msa_file};

  prefetcher_ = std::thread( &read_chunk, std::ref(fptr_), initial_size, std::ref(prefetch_chunk_));
}

const Sequence& MSA_Stream::operator[](const size_t idx) const
{
  if((size_t) idx >= active_chunk_.size())
    throw runtime_error{string("Trying to access MSA_Stream entry out of bounds. i = ") + to_string(idx) };

  return active_chunk_[idx];
}

size_t MSA_Stream::read_next(const size_t number)
{
  // join prefetching thread to ensure new chunk exists
  prefetcher_.join();

  // perform pointer swap to data
  active_chunk_.swap(prefetch_chunk_);

  prefetch_chunk_.clear();
  
  // start request next chunk from prefetcher (async)
  prefetcher_ = std::thread(read_chunk, std::ref(fptr_), number, std::ref(prefetch_chunk_));

  // return size of current buffer
  return active_chunk_.size();
}
