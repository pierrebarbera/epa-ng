#include "MSA_Stream.hpp"

static void read_chunk(MSA_Stream::file_type::pointer fptr, const size_t number, MSA_Stream::container_type& prefetch_buffer)
{
  if (!fptr) {
    throw std::runtime_error{"fptr was invalid!"};
  }

  int sites = 0;
  size_t number_left = number;

  char * sequence = nullptr;
  char * header = nullptr;
  long sequence_length;
  long header_length;
  long sequence_number;

  prefetch_buffer.clear();

  while (number_left and
    pll_fasta_getnext(fptr, &header, &header_length, &sequence, &sequence_length, &sequence_number))
  {
    if (sites && (sites != sequence_length))
      throw std::runtime_error{"MSA file does not contain equal size sequences"};

    if (!sites) sites = sequence_length;

    for (long i = 0; i < sequence_length; ++i) sequence[i] = toupper(sequence[i]);
    
    prefetch_buffer.emplace_back(header, sequence);
    free(sequence);
    free(header);
    sites = sequence_length;
    number_left--;
  }
  printf("read_chunk! returned size: %lu\n", prefetch_buffer.size());
}

MSA_Stream::MSA_Stream (const std::string& msa_file, const size_t initial_size)
  : fptr_(nullptr, fasta_close), active_chunk_(), prefetch_chunk_()
{
  fptr_ = file_type(pll_fasta_open(msa_file.c_str(), pll_map_fasta),
                    fasta_close);
  if(!fptr_) {
    throw std::runtime_error{std::string("Cannot open file: ") + msa_file};
  }

  read_chunk(fptr_.get(), initial_size, prefetch_chunk_);

  // prefetcher_ = std::thread(read_chunk, fptr_.get(), initial_size, std::ref(prefetch_chunk_));
}

const Sequence& MSA_Stream::operator[](const size_t idx) const
{
  if(idx >= active_chunk_.size()) {
    throw std::runtime_error{std::string("Trying to access MSA_Stream entry out of bounds. i = ") + std::to_string(idx) };
  }

  return active_chunk_[idx];
}

size_t MSA_Stream::read_next(const size_t number)
{
  // join prefetching thread to ensure new chunk exists
  if (prefetcher_.joinable()) {
    prefetcher_.join();
  }

  // perform pointer swap to data
  std::swap(active_chunk_, prefetch_chunk_);
  // active_chunk_.swap(prefetch_chunk_);

  if (!fptr_) {
    throw std::runtime_error{"fptr was invalid during read_next!"};
  }
  
  // start request next chunk from prefetcher (async)
  prefetcher_ = std::thread(read_chunk, fptr_.get(), number, std::ref(prefetch_chunk_));

  // read_chunk(fptr_.get(), number, prefetch_chunk_);

  // return size of current buffer
  printf("read_next! returned size: %lu\n", active_chunk_.size());
  return active_chunk_.size();
}

MSA_Stream::~MSA_Stream() 
{
  // avoid dangling threads
  if (prefetcher_.joinable()) {
    prefetcher_.join();
  }
}
