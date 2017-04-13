#include "jplace_util.hpp"

#include <sstream>

using namespace std;

void merge_into(ofstream& dest, const vector<string>& sources)
{
  size_t i = 0;
  for (const auto& file_n : sources)
  {
    ifstream file(file_n);
    dest << file.rdbuf();
    dest.clear(); // empty input files silently set failure flags!
    if (++i < sources.size()) {
      dest << ",";
    }
    dest << NEWL;
    file.close();
  }
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

string pquery_to_jplace_string(const PQuery& pquery, const MSA& msa)
{
  ostringstream output;

  output << "    {\"p\":" << NEWL; // p for pquery
  output << "      [" << NEWL; // opening bracket for pquery array

  size_t i = 0;
  for (const auto& place : pquery)
  {
    // individual pquery
    output << "      " << placement_to_jplace_string(place);
    if (++i < pquery.size()) {
      output << ",";  
    }
    output << NEWL;
  } 

  // closing bracket for pquery array, and start of name column
  output << "      ]," << NEWL <<"    \"n\": [";
  // list of sequence headers
  i = 0;
  for (const auto& header : msa[pquery.sequence_id()].header_list() ) {
    output << "\"" << header.c_str() << "\"";
    if (++i < msa[pquery.sequence_id()].header_list().size()) {
      output << ",";  
    }
  }

  output << "]" << NEWL; // close name bracket

  output << "    }";// final bracket

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

string sample_to_jplace_string(const Sample& sample, const MSA& msa)
{
  ostringstream output;

  size_t i = 0;
  for (const auto& p : sample) {
    output << pquery_to_jplace_string(p, msa); 
    if (++i < sample.size()) {
      output << ",";
    }
    output << NEWL;
  }
  return output.str();
}

string full_jplace_string(const Sample& sample, 
                          const string& invocation, 
                          const MSA& msa)
{
  ostringstream output;

  // tree and other init
  output << init_jplace_string(sample.newick());

  // actual placements
  output << sample_to_jplace_string(sample, msa);

  // metadata string
  output << finalize_jplace_string(invocation);

  return output.str();
}
