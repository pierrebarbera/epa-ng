#include "jplace_util.hpp"

using namespace std;

void merge_into(ofstream& dest, vector<string>& sources)
{
  for(auto& file_n : sources)
  {
    ifstream file(file_n);
    dest << file.rdbuf();
    rwnd(dest, 1);
    dest << "," << NEWL;
    file.close();
  }
  rwnd(dest, 2);
  dest << NEWL;
}

string placement_to_jplace_string(const Placement& p)
{
  ostringstream output;

  output << "[" << to_string(p.branch_id()) << ", ";
  output << to_string(p.likelihood()) << ", ";
  output << to_string(p.lwr()) << ", ";
  output << to_string(p.distal_length()) << ", ";
  output << to_string(p.pendant_length()) << "]";

  return output.str();
}

string init_jplace_string(const string& numbered_newick)
{
  ostringstream output;

  output << "{" << NEWL;
  output << "  \"tree\": \"" << numbered_newick << "\"," << NEWL;
  output << "  \"placements\": " << NEWL;
  output << "  [" << NEWL;

  return output.str();
}

string finalize_jplace_string(const string& invocation)
{
  assert(invocation.length() > 0);

  ostringstream output;

  output << "  ]," << NEWL;

  output << "  \"metadata\": {\"invocation\": \"" << invocation << "\"}," << NEWL;

  output << "  \"version\": 3," << NEWL;
  output << "  \"fields\": ";
  output << "[\"edge_num\", \"likelihood\", \"like_weight_ratio\", \"distal_length\"";
  output << ", \"pendant_length\"]" << NEWL;

  output << "}" << NEWL;

  return output.str();
}
