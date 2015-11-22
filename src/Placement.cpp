#include "Placement.hpp"

#include <stdexcept>

using namespace std;

Placement::Placement ()
  : likelihoods_(0), sequence(NULL)
{ }

Placement::Placement (const unsigned int size, const Sequence * s)
  : likelihoods_(size), sequence(s)
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
