#include "src/epa.hpp"
#include "Epatest.hpp"

TEST(epa, epa)
{
  std::string tree_file(env->data_dir);
  tree_file += "ref.tre";
  std::string reference_file(env->data_dir);
  reference_file += "aln.fasta";
  std::string query_file(env->data_dir);
  query_file += "query.fasta";

  epa(tree_file,
      reference_file,
      query_file,
      {0.25, 0.25, 0.25, 0.25},
      {1,1,1,1,1,1});
}
