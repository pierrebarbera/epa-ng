#include "Epatest.hpp"

#include "src/pllhead.hpp"
#include "src/file_io.hpp"
#include "src/epa_pll_util.hpp"
#include "src/jplace_util.hpp"
#include "src/Tree.hpp"
#include "src/MSA.hpp"

#include <string>
#include <vector>
#include <iostream>

using namespace std;

TEST(jplace_util, placement_to_jplace_string)
{
  // buildup
  auto tree = Tree(env->tree_file, env->reference_file, env->model);
  auto query_reads = build_MSA_from_file(env->query_file);

  // tests
  auto placements = tree.place(query_reads);
  vector<string> out;

  for (auto const &p : placements)
    out.push_back(placement_to_jplace_string(p));

  EXPECT_EQ(out.size(), 2);

  // teardown

}
