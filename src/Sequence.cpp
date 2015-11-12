#include "Sequence.hpp"

#include <string>

using namespace std;

Sequence::Sequence(const Sequence& s)
  : header_(s.header())
  , sequence_(s.sequence())
{}

Sequence::Sequence(string header, string sequence)
          : header_(header), sequence_(sequence)
{

}

Sequence::~Sequence()
{

}
