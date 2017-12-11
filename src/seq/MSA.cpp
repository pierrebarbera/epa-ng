#include "seq/MSA.hpp"

#include <stdexcept>

void MSA::move_sequences(MSA::iterator begin, MSA::iterator end)
{
  std::move(begin, end, std::back_inserter(sequence_list_));
}

const Sequence& MSA::operator[](const size_t i) const
{
  //if((size_t) i >= sequence_list_.size()) {
  //  throw std::runtime_error{std::string("Trying to access MSA entry out of bounds. i = ") + std::to_string(i) };
  //}

  return sequence_list_[i];
}

void MSA::append(const std::string& header, const std::string& sequence)
{
  if(num_sites_ && sequence.length() != num_sites_) {
    throw std::runtime_error{std::string("Tried to insert sequence to MSA of unequal length: ") + header};
  }

  sequence_list_.emplace_back(header, sequence);

  if (!num_sites_) {
    num_sites_ = sequence.length();
  }
}

void std::swap(MSA& a, MSA& b)
{
  MSA::swap(a, b);
}
