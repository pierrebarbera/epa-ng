#include "jplace_util.hpp"

using namespace std;

static void rwnd(ostringstream& ss, unsigned int chars)
{
  long pos = ss.tellp();
  ss.seekp (pos - (chars * sizeof(char)));
}

string pquery_to_jplace_string(const PQuery& p)
{
  ostringstream output;

  output << "    {\"p\":" << NEWL; // p for pquery
  output << "      [" << NEWL; // opening bracket for pquery array

  unsigned int branch = 0; // TODO bad, will change when map introduced
  for (auto logl : p)
  {
    // individual pquery
    output << "      [" << branch << ", "<< logl << "]," << NEWL;
    ++branch;
  }

  // rewind by two chars: newline and last comma
  rwnd(output, 2);

  // closing bracket for pquery array and name column
  output << "\n      ],\n    \"n\": [\"" << p.sequence().header().c_str() << "\"]" << NEWL;

  output << "    }";// final bracket

  return output.str();
}

std::string pquery_set_to_jplace_string(const PQuery_Set& ps, string& invocation)
{
  ostringstream output;

  output << "{\n";

  output << "  \"tree\": \"" << ps.newick() << "\",\n";
  output << "  \"pquerys\": \n";
  output << "  [\n";

  for (auto p : ps)
    output << pquery_to_jplace_string(p) << ",\n";

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
