#include "Epatest.hpp"

#include "src/file_io.hpp"
#include "src/Tree.hpp"
#include "src/MSA.hpp"
#include "src/place.hpp"

#include <string>
#include <vector>
#include <limits>

using namespace std;

TEST(Tree, process_from_binary)
{
  #ifndef __MPI
  // setup
  auto msa = build_MSA_from_file(env->reference_file);
  auto queries = MSA_Stream(env->query_file);
  Model model;
  Options options;
  Tree original_tree(env->tree_file, msa, model, options);
  dump_to_binary(original_tree, env->binary_file);
  Tree read_tree(env->binary_file, model, options);

  EXPECT_DOUBLE_EQ(original_tree.ref_tree_logl(), read_tree.ref_tree_logl());

  // test
  string invocation("./this --is -a test");

  process(read_tree, queries, env->out_dir, options, invocation);

  options.prescoring = true;
  process(read_tree, queries, env->out_dir, options, invocation);

  Tree mvstree;
  mvstree = Tree(env->binary_file, model, options);

  process(mvstree, queries, env->out_dir, options, invocation);

  // teardown
  #endif
}

TEST(Tree, combined_input_file)
{
  auto combined_msa = build_MSA_from_file(env->combined_file);
  auto tree = Tree(env->tree_file, combined_msa, env->model, env->options);
}

TEST(Tree, rooted_input)
{
  auto msa = build_MSA_from_file(env->reference_file);
  auto tree = Tree(env->tree_file_rooted, msa, env->model, env->options);
}
