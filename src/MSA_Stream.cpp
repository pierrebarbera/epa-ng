#include "MSA_Stream.hpp"

using namespace std;

MSA_Stream::MSA_Stream (const string& msa_file)
{
  fptr_ = pll_fasta_open(msa_file.c_str(), pll_map_fasta);
  if(!fptr_)
    throw runtime_error{string("Cannot open file: ") + msa_file};
}

const Sequence& MSA_Stream::operator[](const unsigned int idx) const
{
  assert(fptr_);
  if((size_t) idx >= sequence_list_.size())
    throw runtime_error{string("Trying to access MSA entry out of bounds. i = ") + to_string(idx) };

  return sequence_list_[idx];
}

unsigned int MSA_Stream::read_next(const unsigned int number)
{
  assert(fptr_);
  int sites = 0;
  unsigned int number_left = number;

  char * sequence = NULL;
  char * header = NULL;
  long sequence_length;
  long header_length;
  long sequence_number;

  while (number_left and pll_fasta_getnext(fptr_, &header, &header_length, &sequence, &sequence_length, &sequence_number))
  {
    if (sites && (sites != sequence_length))
      throw runtime_error{"MSA file does not contain equal size sequences"};

    if (!sites) sites = sequence_length;

    sequence_list_.emplace_back(header, sequence);
    free(sequence);
    free(header);
    sites = sequence_length;
    number_left--;
  }
  return number - number_left;
}
