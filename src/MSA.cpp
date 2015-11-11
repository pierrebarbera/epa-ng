#include "MSA.h"

#include <vector>
#include <tuple>
#include <stdexcept>

using namespace std;

// Invariant: all contained sequences shall have the same length (num_sites)
MSA::MSA(const int num_sites) : num_sites(num_sites)
{
  if(num_sites <= 0)
    throw runtime_error{"Number of sites in MSA cannot be zero or negative."};
}

MSA::~MSA()
{
  // vectors are self deleting
}


tuple<string, string> MSA::get(const int i) const
{
  if(i >= sequence_list.size() || i < 0)
    throw runtime_error{"Trying to access MSA entry out of bounds"};

  return make_tuple(sequence_list[i].header, sequence_list[i].sequence);
}

template <typename T, typename D>
void MSA::append(T header, D sequence)
{
  if(sequence.length() != num_sites)
    throw runtime_error{string("Tried to insert sequence to MSA of unequal length: ") + sequence};

  sequence_list.push_back(Sequence(header, sequence));
}

MSA::iterator MSA::begin()
{
    return sequence_list.begin();
}

MSA::iterator MSA::end()
{
    return sequence_list.end();
}

MSA::const_iterator MSA::cbegin()
{
    return sequence_list.cbegin();
}

MSA::const_iterator MSA::cend()
{
    return sequence_list.cend();
}
