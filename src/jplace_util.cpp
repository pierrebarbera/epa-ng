#include "jplace_util.hpp"

#include <sstream>

using namespace std;

string placement_to_jplace_string(const Placement& p)
{
  ostringstream output;

  output << "    {\"p\":" << NEWL; // p for placement
  output << "      [" << NEWL; // opening bracket for placement array

  unsigned int branch = 0; // TODO bad, will change when map introduced
  for (auto logl : p)
  {
    // individual placement
    output << "      [" << branch << ", "<< logl << "]," << NEWL;
    ++branch;
  }

  // rewind by two chars: newline and last comma
  long pos = output.tellp();
  output.seekp (pos - (2 * sizeof(char)));

  // closing bracket for placement array and name column
  output << "\n      ],\n    \"n\": [\"" << p.sequence().header().c_str() << "\"]" << NEWL;

  output << "    }" << NEWL;// final bracket

  return output.str();
}

std::string placement_set_to_jplace_string(const Placement_Set& ps)
{
  ostringstream output;

  output << "{\n";

  // TODO newick string with edge numbering
  output << "  \"tree\": \"" << "NEWICK-NEWICK" << "\",\n";
  output << "  \"placements\": \n";
  output << "  [\n";

  for (auto p : ps)
    output << placement_to_jplace_string(p);

  output << "  ],\n";

  // TODO pass metadata string to this function?

  output << "  \"version\": 3,\n";
  output << "  \"fields\": \n";
  output << "  [\"edge_num\", \"likelihood\"]\n";

  output << "}\n";

  return output.str();
}
