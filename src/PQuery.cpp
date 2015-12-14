#include "PQuery.hpp"

#include <stdexcept>

using namespace std;

PQuery::PQuery (const Sequence& s)
  : sequence_(s)
{ }

PQuery::~PQuery ()
{ }

PQuery::iterator PQuery::begin()
{
    return placements_.begin();
}

PQuery::iterator PQuery::end()
{
    return placements_.end();
}


PQuery::const_iterator PQuery::begin() const
{
    return placements_.cbegin();
}

PQuery::const_iterator PQuery::end() const
{
    return placements_.cend();
}

PQuery::const_iterator PQuery::cbegin()
{
    return placements_.cbegin();
}

PQuery::const_iterator PQuery::cend()
{
    return placements_.cend();
}
