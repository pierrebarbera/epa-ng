#include "MSA.hpp"

#include <vector>
#include <tuple>
#include <stdexcept>

#include "Sequence.hpp"

using namespace std;

void MSA::move_sequences(iterator begin, iterator end)
{
  // sequence_list_.assign(begin, end);
  std::move(begin, end, std::back_inserter(sequence_list_));
}

const Sequence& MSA::operator[](const unsigned int i) const
{
  if((size_t) i >= sequence_list_.size())
    throw runtime_error{string("Trying to access MSA entry out of bounds. i = ") + to_string(i) };

  return sequence_list_[i];
}

void MSA::append(const string& header, const string& sequence)
{
  if(num_sites_ && sequence.length() != (size_t) num_sites_)
    throw runtime_error{string("Tried to insert sequence to MSA of unequal length: ") + sequence};

  sequence_list_.emplace_back(header, sequence);

  if (!num_sites_)
    num_sites_ = sequence.length();
}

void MSA::append(Sequence s)
{
  if(num_sites_ && s.sequence().length() != (size_t) num_sites_)
    throw runtime_error{string("Tried to insert sequence to MSA of unequal length: ") + s.sequence()};

  sequence_list_.emplace_back(s);

  if (!num_sites_)
    num_sites_ = s.sequence().length();

}
