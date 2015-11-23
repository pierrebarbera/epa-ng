#include "jplace_util.hpp"

using namespace std;

static void rwnd(ostringstream& ss, unsigned int chars)
{
  long pos = ss.tellp();
  ss.seekp (pos - (chars * sizeof(char)));
}

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
  rwnd(output, 2);

  // closing bracket for placement array and name column
  output << "\n      ],\n    \"n\": [\"" << p.sequence().header().c_str() << "\"]" << NEWL;

  output << "    }";// final bracket

  return output.str();
}

std::string placement_set_to_jplace_string(const Placement_Set& ps, string& invocation)
{
  ostringstream output;

  output << "{\n";

  output << "  \"tree\": \"" << ps.newick() << "\",\n";
  output << "  \"placements\": \n";
  output << "  [\n";

  for (auto p : ps)
    output << placement_to_jplace_string(p) << ",\n";

  // undo the last comma
  rwnd(output, 2);

  output << "\n  ],\n";

  // metadata string
  if (invocation.length() > 0)
    output << "\"metadata\": {\"invocation\": \"" << invocation << "\"},\n";

  output << "  \"version\": 3,\n";
  output << "  \"fields\": ";
  output << "[\"edge_num\", \"likelihood\"]\n";

  output << "}\n";

  return output.str();
}
