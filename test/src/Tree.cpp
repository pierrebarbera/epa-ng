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

  EXPECT_DOUBLE_EQ(original_tree.ref_tree_logl(), read_tree.ref_tree_logl());

  // test
  string invocation("./this --is -a test");

  process(read_tree, queries, env->out_dir, options, invocation);

  options.prescoring = true;
  process(read_tree, queries, env->out_dir, options, invocation);

  Tree mvstree;
  mvstree = Tree(env->binary_file, options);

  process(mvstree, queries, env->out_dir, options, invocation);

  // teardown
}

TEST(Tree, place_from_binary)
{
  // setup
  auto msa = build_MSA_from_file(env->reference_file);
  auto queries = build_MSA_from_file(env->query_file);
  Model model;
  Options options;
  options.support_threshold = 1.0;
  options.acc_threshold = true;
  Tree original_tree(env->tree_file, msa, model, options);
  dump_to_binary(original_tree, env->binary_file);
  Tree read_tree(env->binary_file, options);


  // test
  ASSERT_DOUBLE_EQ(original_tree.options().support_threshold, read_tree.options().support_threshold);
  ASSERT_DOUBLE_EQ(original_tree.ref_tree_logl(), read_tree.ref_tree_logl());

  auto original_sample = place(original_tree, queries);
  auto read_sample = place(read_tree, queries);

  EXPECT_DOUBLE_EQ(original_tree.ref_tree_logl(), read_tree.ref_tree_logl());

  ASSERT_EQ(original_sample.size(), read_sample.size());

  for (size_t i = 0; i < original_sample.size(); i++)
  {
    auto& o_pq = original_sample[i];

    // find the correct pquery
    PQuery r_pq;
    bool found_pq = false;
    for(auto& pq : read_sample)
    {
      if (pq.sequence_id() == o_pq.sequence_id())
      {
        found_pq = true;
        r_pq = pq;
      }
    }
    ASSERT_TRUE(found_pq);

    ASSERT_EQ(o_pq.size(), r_pq.size());

    for (size_t j = 0; j < o_pq.size(); j++)
    {
      auto original_place = o_pq[j];
      // find correct placement
      Placement read_place;
      bool found_p = false;
      for(auto& p : r_pq)
      {
        if (p.branch_id() == original_place.branch_id())
        {
          found_p = true;
          read_place = p;
        }
      }
      ASSERT_TRUE(found_p);

      EXPECT_DOUBLE_EQ(original_place.likelihood(), read_place.likelihood());
      // printf("edgenum %d/%d: %f vs %f\n",original_place.branch_id(), read_place.branch_id(),
      // original_place.likelihood(), read_place.likelihood());
    }
  }
  // teardown
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
