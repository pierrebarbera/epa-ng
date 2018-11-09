#include "io/jplace_util.hpp"

#include <sstream>
#include <tuple>

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

void placement_to_jplace_string(Placement const& p,
                                std::ostream& os,
                                rtree_mapper const& mapper)
{
  auto branch_id = p.branch_id();
  auto distal_length = p.distal_length();
  if ( mapper ) {
    std::tie(branch_id, distal_length) = mapper.in_rtree(branch_id, distal_length);
  }

  os << "[" << branch_id << ", ";
  os << p.likelihood() << ", ";
  os << p.lwr() << ", ";
  os << distal_length << ", ";
  os << p.pendant_length() << "]";
}

void pquery_to_jplace_string( PQuery<Placement> const& pquery,
                              std::ostream& os,
                              rtree_mapper const& mapper)
{
  os << "    {\"p\": [" << NEWL; // p for pquery

  size_t i = 0;
  for (const auto& place : pquery)
  {
    // individual pquery
    os << "      ";

    placement_to_jplace_string(place, os, mapper);

    if (++i < pquery.size()) {
      os << ",";
    }
    os << NEWL;
  }

  // closing bracket for pquery array
  os << "      ]," << NEWL;

  // start of name column
  os <<"    \"n\": [";

  // sequence header
  const auto& header = pquery.header();
  os << "\"" << header.c_str() << "\"";


  os << "]" << NEWL; // close name bracket

  os << "    }";// final bracket

}

void init_jplace_string(const std::string& numbered_newick, std::ostream& os)
{
  os << "{" << NEWL;
  os << "  \"tree\": \"" << numbered_newick << "\"," << NEWL;
  os << "  \"placements\": " << NEWL;
  os << "  [" << NEWL;
}

void finalize_jplace_string(const std::string& invocation, std::ostream& os)
{
  assert(invocation.length() > 0);

  os << "  ]," << NEWL;

  os << "  \"metadata\": {\"invocation\": \"" << invocation << "\"}," << NEWL;

  os << "  \"version\": 3," << NEWL;
  os << "  \"fields\": ";
  os << "[\"edge_num\", \"likelihood\", \"like_weight_ratio\", \"distal_length\"";
  os << ", \"pendant_length\"]" << NEWL;

  os << "}" << NEWL;
}

void sample_to_jplace_string(Sample<Placement> const& sample, std::ostream& os, rtree_mapper const& mapper)
{
  size_t i = 0;
  for (const auto& p : sample) {
    pquery_to_jplace_string(p, os, mapper);
    if (++i < sample.size()) {
      os << ",";
    }
    os << NEWL;
  }
}

void full_jplace_string(Sample<Placement> const& sample,
                        std::string const& invocation,
                        std::ostream& os,
                        rtree_mapper const& mapper)
{
  // tree and other init
  init_jplace_string(sample.newick(), os);

  // actual placements
  sample_to_jplace_string(sample, os, mapper);

  // metadata std::string
  finalize_jplace_string(invocation, os);

}
