#include "jplace_util.hpp"

#include <sstream>

using namespace std;

string placement_to_jplace_string(const Placement& p)
{
  ostringstream output;

  output << "[" << to_string(p.branch_id()) << ", ";
  output << to_string(p.likelihood()) << ", ";
  output << to_string(p.lwr()) << ", ";
  output << to_string(p.pendant_length()) << ", ";
  output << to_string(p.distal_length()) << "]";

  return output.str();
}

string pquery_to_jplace_string(const PQuery& pquery)
{
  ostringstream output;

  output << "    {\"p\":" << NEWL; // p for pquery
  output << "      [" << NEWL; // opening bracket for pquery array

  for (auto place : pquery)
  {
    // individual pquery
    output << "      " << placement_to_jplace_string(place) << "," << NEWL;
  }

  // undo last comma and newline
  rwnd(output, 2);
  output << NEWL;

  // closing bracket for pquery array, and start of name column
  output << "      ]," << NEWL <<"    \"n\": [";
  // list of sequence headers
  for (const auto& header : pquery.sequence().header_list() )
    output << "\"" << header.c_str() << "\",";

  // rewind last comma
  rwnd(output, 1);

  output << "]" << NEWL; // close name bracket

  output << "    }";// final bracket

  return output.str();
}

std::string pquery_set_to_jplace_string(const PQuery_Set& ps, string& invocation)
{
  ostringstream output;

  output << "{" << NEWL;

  output << "  \"tree\": \"" << ps.newick() << "\"," << NEWL;
  output << "  \"placements\": " << NEWL;
  output << "  [" << NEWL;

  for (auto p : ps)
    output << pquery_to_jplace_string(p) << "," << NEWL;

  // undo the last comma
  rwnd(output, 2);
  output << NEWL;

  output << "  ]," << NEWL;

  // metadata string
  if (invocation.length() > 0)
    output << "  \"metadata\": {\"invocation\": \"" << invocation << "\"}," << NEWL;

  output << "  \"version\": 3," << NEWL;
  output << "  \"fields\": ";
  output << "[\"edge_num\", \"likelihood\", \"like_weight_ratio\", \"pendant_length\"";
  output << ", \"distal_length\"]" << NEWL;

  output << "}" << NEWL;

  return output.str();
}
