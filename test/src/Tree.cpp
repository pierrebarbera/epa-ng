#include "Epatest.hpp"

#include "src/file_io.hpp"
#include "src/Tree.hpp"
#include "src/MSA.hpp"
#include "src/place.hpp"

#include <string>
#include <vector>

using namespace std;

TEST(Tree, place)
{
  // buildup
  auto query_msa = build_MSA_from_file(env->query_file);
  auto reference_msa = build_MSA_from_file(env->reference_file);
  Options options;
  options.support_threshold = 0.0;
  auto tree = Tree(env->tree_file, reference_msa, env->model, options);

  // tests
  auto sample = place(tree, query_msa);

  EXPECT_EQ(sample.size(), 2);

  // check for zero entires and correct number
  for (auto pquery : sample)
  {
    int count = 0;
    for (auto p : pquery)
    {
      EXPECT_NE(p.likelihood(), 0.0);
      count++;
    }
    EXPECT_EQ(count, tree.nums().branches);
  }
  // teardown

}

TEST(Tree, place_prescore)
{
  // buildup
  auto query_msa = build_MSA_from_file(env->query_file);
  auto reference_msa = build_MSA_from_file(env->reference_file);
  Options options;
  options.prescoring = true;
  options.opt_model = true;
  options.opt_branches = true;
  options.support_threshold = 0.0;
  auto tree = Tree(env->tree_file, reference_msa, env->model, options);

  // tests
  auto sample = place(tree, query_msa);

  EXPECT_EQ(sample.size(), 2);

  // check for zero entires and correct number
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
  // teardown

}

TEST(Tree, process_from_binary)
{
  // setup
  auto msa = build_MSA_from_file(env->reference_file);
  auto queries = MSA_Stream(env->query_file);
  Model model;
  Options options;
  Tree original_tree(env->tree_file, msa, model, options);
  dump_to_binary(original_tree, env->binary_file);
  Tree read_tree(env->binary_file, options);

  // test
  string invocation("./this --is -a test");

  process(read_tree, queries, env->out_dir, options, invocation);

  // teardown
}

TEST(Tree, combined_input_file)
{
  auto combined_msa = build_MSA_from_file(env->combined_file);
  auto tree = Tree(env->tree_file, combined_msa, env->model, env->options);
}
