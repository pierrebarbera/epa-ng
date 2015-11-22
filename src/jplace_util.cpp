#include "jplace_util.hpp"

#include <sstream>

using namespace std;

string placement_to_jplace_string(const Placement& p)
{
  ostringstream output;

  output << "{\"p\":" << NEWL; // p for placement
  output << "  [" << NEWL; // opening bracket for placement array

  unsigned int branch = 0; // TODO bad, will change when map introduced
  for (auto logl : p)
  {
    // individual placement
    output << "  [" << branch << ", "<< logl << "]," << NEWL;
    ++branch;
  }

  // rewind by two chars: newline and last comma
  long pos = output.tellp();
  output.seekp (pos - (2 * sizeof(char)));

  // closing bracket for placement array and name column
  output << "\n  ],\n\"n\": [\"" << p.sequence().header().c_str() << "\"]" << NEWL;

  output << "}" << NEWL;// final bracket

  return output.str();
}
