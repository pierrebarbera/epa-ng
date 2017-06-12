#include "jplace_util.hpp"

#include <sstream>

void merge_into(std::ofstream& dest, const std::vector<std::string>& sources)
{
  size_t i = 0;
  for (const auto& file_n : sources)
  {
    std::ifstream file(file_n);
    dest << file.rdbuf();
    dest.clear(); // empty input files silently set failure flags!
    if (++i < sources.size()) {
      dest << ",";
    }
    dest << NEWL;
    file.close();
  }
}

std::string placement_to_jplace_string(const Placement& p)
{
  std::ostringstream output;

  output << "[" << std::to_string(p.branch_id()) << ", ";
  output << std::to_string(p.likelihood()) << ", ";
  output << std::to_string(p.lwr()) << ", ";
  output << std::to_string(p.distal_length()) << ", ";
  output << std::to_string(p.pendant_length()) << "]";

  return output.str();
}

std::string pquery_to_jplace_string(const PQuery& pquery, const MSA& msa)
{
  std::ostringstream output;

  output << "    {\"p\": [" << NEWL; // p for pquery

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

  // closing bracket for pquery array
  output << "      ]," << NEWL; 
  
  // start of entropy column
  output << "    \"m\": {\"entropy\": " << std::to_string(pquery.entropy()) << "}," << NEWL;

  // start of name column
  output <<"    \"n\": [";
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

std::string init_jplace_string(const std::string& numbered_newick)
{
  std::ostringstream output;

  output << "{" << NEWL;
  output << "  \"tree\": \"" << numbered_newick << "\"," << NEWL;
  output << "  \"placements\": " << NEWL;
  output << "  [" << NEWL;

  return output.str();
}

std::string finalize_jplace_string(const std::string& invocation)
{
  assert(invocation.length() > 0);

  std::ostringstream output;

  output << "  ]," << NEWL;

  output << "  \"metadata\": {\"invocation\": \"" << invocation << "\"}," << NEWL;

  output << "  \"version\": 3," << NEWL;
  output << "  \"fields\": ";
  output << "[\"edge_num\", \"likelihood\", \"like_weight_ratio\", \"distal_length\"";
  output << ", \"pendant_length\"]" << NEWL;

  output << "}" << NEWL;

  return output.str();
}

std::string sample_to_jplace_string(const Sample& sample, const MSA& msa)
{
  std::ostringstream output;

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

std::string full_jplace_string(const Sample& sample, 
                          const std::string& invocation, 
                          const MSA& msa)
{
  std::ostringstream output;

  // tree and other init
  output << init_jplace_string(sample.newick());

  // actual placements
  output << sample_to_jplace_string(sample, msa);

  // metadata std::string
  output << finalize_jplace_string(invocation);

  return output.str();
}
