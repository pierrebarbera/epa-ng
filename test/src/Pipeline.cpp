#include "Epatest.hpp"

#include "src/Pipeline_Stage.hpp"
#include "src/MSA.hpp"
#include "src/MSA_Stream.hpp"
#include "src/Tree.hpp"
#include "src/Model.hpp"
#include "src/Options.hpp"
#include "src/pll_util.hpp"
#include "src/file_io.hpp"
#include "src/Sample.hpp"

using namespace std;

TEST(Pipeline, memory)
{
  auto msa = build_MSA_from_file(env->reference_file);
  auto queries = MSA_Stream(env->query_file);
  Model model;
  Options options;
  Tree tree(env->tree_file, msa, model, options);

  // output class
  Sample sample;

  // Prepare tiny trees
  const auto num_branches = tree.nums().branches;
  vector<pll_utree_t *> branches(num_branches);
  auto num_traversed_branches = utree_query_branches(tree.tree(), &branches[0]);
  assert(num_traversed_branches == num_branches);

  // build all tiny trees with corresponding edges
  vector<Tiny_Tree> insertion_trees;
  for (unsigned int branch_id = 0; branch_id < num_branches; branch_id++)
  {
    insertion_trees.emplace_back(branches[branch_id], branch_id, tree, !options.prescoring);
  }

  // set up the pipeline
  Memory_Communicator<MSA_Stream*> ingestion(&queries);
  Memory_Communicator<Sample*> output(&sample);

  Placement_Stage place_stage(&ingestion, &output, &insertion_trees);

  // run
  place_stage.process();

  // some tests
  for (auto pquery : sample)
  {
    int count = 0;
    for (auto p : pquery)
    {
      EXPECT_NE(p.likelihood(), 0.0);
      count++;
    }
    EXPECT_EQ(count, 1);
  }

}
