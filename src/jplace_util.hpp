#pragma once

#include <fstream>
#include <vector>
#include <sstream>

#include "stringify.hpp"
#include "PQuery.hpp"
#include "Sample.hpp"
#include "Placement.hpp"
#include "MSA_Stream.hpp"

std::string placement_to_jplace_string(const Placement& p);
std::string init_jplace_string(const std::string& numbered_newick);
std::string finalize_jplace_string(const std::string& invocation);
void merge_into(std::ofstream& dest, std::vector<std::string>& sources);

template<typename msa_type>
std::string pquery_to_jplace_string(const PQuery& pquery, const msa_type& msa)
{
  std::ostringstream output;

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
  for (const auto& header : msa[pquery.sequence_id()].header_list() )
    output << "\"" << header.c_str() << "\",";

  // rewind last comma
  rwnd(output, 1);

  output << "]" << NEWL; // close name bracket

  output << "    }";// final bracket

  return output.str();
}

template<typename msa_type>
std::string sample_to_jplace_string(const Sample& sample, const msa_type& msa)
{
  std::ostringstream output;

  for (auto p : sample)
    output << pquery_to_jplace_string(p, msa) << "," << NEWL;

  // undo the last comma
  rwnd(output, 2);
  output << NEWL;
  return output.str();
}

template<typename msa_type>
std::string full_jplace_string(const Sample& sample, const std::string& invocation, const msa_type& msa)
{
  std::ostringstream output;

  // tree and other init
  output << init_jplace_string(sample.newick());

  // actual placements
  output << sample_to_jplace_string(sample, msa);

  // metadata string
  output << finalize_jplace_string(invocation);

  return output.str();
}
