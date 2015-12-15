#include "Epatest.hpp"

#include "src/file_io.hpp"
#include "src/Tree.hpp"
#include "src/MSA.hpp"

#include <string>
#include <vector>

using namespace std;

TEST(Tree, place)
{
  // buildup
  auto query_msa = build_MSA_from_file(env->query_file);
  auto reference_msa = build_MSA_from_file(env->reference_file);
  auto tree = Tree(env->tree_file, reference_msa, env->model, env->options, query_msa);

  // tests
  auto pquery_set = tree.place();

  EXPECT_EQ(pquery_set.size(), 2);

  // check for zero entires and correct number
  for (auto pquery : pquery_set)
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

TEST(Tree, place_optimize)
{
  // buildup
  auto query_msa = build_MSA_from_file(env->query_file);
  auto reference_msa = build_MSA_from_file(env->reference_file);
  Options options(true,true,true);
  auto tree = Tree(env->tree_file, reference_msa, env->model, options, query_msa);

  // tests
  auto pquery_set = tree.place();

  EXPECT_EQ(pquery_set.size(), 2);

  // check for zero entires and correct number
  for (auto pquery : pquery_set)
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

TEST(Tree, combined_input_file)
{
  auto combined_msa = build_MSA_from_file(env->combined_file);
  auto tree = Tree(env->tree_file, combined_msa, env->model, env->options);
}
