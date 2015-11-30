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
  auto tree = Tree(env->tree_file, reference_msa, env->model, query_msa);

  // tests
  auto placements = tree.place();

  EXPECT_EQ(placements.size(), 2);

  // check for zero entires and correct number
  for (auto place : placements)
  {
    int count = 0;
    for (auto logl : place)
    {
      EXPECT_NE(logl, 0.0);
      count++;
    }
    EXPECT_EQ(count, tree.nums().branches);
  }
  // teardown

}

TEST(Tree, combined_input_file)
{
  auto combined_msa = build_MSA_from_file(env->combined_file);
  auto tree = Tree(env->tree_file, combined_msa, env->model);
}
