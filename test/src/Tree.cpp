#include "Epatest.hpp"

#include "src/pllhead.hpp"
#include "src/file_io.hpp"
#include "src/jplace_util.hpp"
#include "src/Tree.hpp"
#include "src/MSA.hpp"

#include <string>
#include <vector>

using namespace std;

TEST(Tree, place)
{
  // buildup
  auto tree = Tree(env->tree_file, env->reference_file, env->model);
  auto query_reads = build_MSA_from_file(env->query_file);

  // tests
  auto placements = tree.place(query_reads);

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
