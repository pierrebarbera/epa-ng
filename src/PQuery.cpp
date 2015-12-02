#include "PQuery.hpp"

#include <stdexcept>

using namespace std;

PQuery::PQuery ()
  : likelihoods_(0), sequence_("","")
{ }

PQuery::PQuery (const unsigned int size, const Sequence& s)
  : likelihoods_(size), sequence_(s)
{ }

PQuery::~PQuery ()
{ }

void PQuery::set(const unsigned int branch, const double logl)
{
  if (likelihoods_.size() <= branch)
    throw runtime_error{string("Invalid pquery location: ") + to_string(branch)};
  else
    likelihoods_[branch] = logl;
}

PQuery::iterator PQuery::begin()
{
    return likelihoods_.begin();
}

PQuery::iterator PQuery::end()
{
    return likelihoods_.end();
}


PQuery::const_iterator PQuery::begin() const
{
    return likelihoods_.cbegin();
}

PQuery::const_iterator PQuery::end() const
{
    return likelihoods_.cend();
}

PQuery::const_iterator PQuery::cbegin()
{
    return likelihoods_.cbegin();
}

PQuery::const_iterator PQuery::cend()
{
    return likelihoods_.cend();
}
