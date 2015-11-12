#include "MSA.hpp"

#include <vector>
#include <tuple>
#include <stdexcept>

#include "Sequence.hpp"

using namespace std;

// Invariant: all contained sequences shall have the same length (num_sites)
MSA::MSA(const int num_sites) : num_sites_(num_sites)
{
  if(num_sites <= 0)
    throw runtime_error{"Number of sites in MSA cannot be zero or negative."};
}

MSA::MSA() : num_sites_(0)
{

}

MSA::~MSA()
{
  // vectors are self deleting
}


const Sequence& MSA::get(const int i) const
{
  if((size_t) i >= sequence_list_.size() || i < 0)
    throw runtime_error{string("Trying to access MSA entry out of bounds. i = ") + to_string(i) };

  return sequence_list_[i];
}

void MSA::append(const string& header, const string& sequence)
{
  if(sequence.length() != (size_t) num_sites_)
    throw runtime_error{string("Tried to insert sequence to MSA of unequal length: ") + sequence};

  // sequence_list.push_back(Sequence(header, sequence));
  sequence_list_.emplace_back(header, sequence);
}

void MSA::append(Sequence s)
{
  sequence_list_.push_back(s);
}

MSA::iterator MSA::begin()
{
    return sequence_list_.begin();
}

MSA::iterator MSA::end()
{
    return sequence_list_.end();
}

MSA::const_iterator MSA::cbegin()
{
    return sequence_list_.cbegin();
}

MSA::const_iterator MSA::cend()
{
    return sequence_list_.cend();
}
