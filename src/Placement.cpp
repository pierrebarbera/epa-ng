#include "Placement.hpp"

#include <stdexcept>

using namespace std;

Placement::Placement ()
  : likelihoods_(0), sequence_("","")
{ }

Placement::Placement (const unsigned int size, const Sequence& s)
  : likelihoods_(size), sequence_(s)
{ }

Placement::~Placement ()
{ }

void Placement::set(const unsigned int branch, const double logl)
{
  if (likelihoods_.size() <= branch)
    throw runtime_error{string("Invalid placement location: ") + to_string(branch)};
  else
    likelihoods_[branch] = logl;
}

Placement::iterator Placement::begin()
{
    return likelihoods_.begin();
}

Placement::iterator Placement::end()
{
    return likelihoods_.end();
}


Placement::const_iterator Placement::begin() const
{
    return likelihoods_.cbegin();
}

Placement::const_iterator Placement::end() const
{
    return likelihoods_.cend();
}

Placement::const_iterator Placement::cbegin()
{
    return likelihoods_.cbegin();
}

Placement::const_iterator Placement::cend()
{
    return likelihoods_.cend();
}
