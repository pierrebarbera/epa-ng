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

TEST(jplace_util, pquery_to_jplace_string)
{
  // buildup
  auto query_msa = build_MSA_from_file(env->query_file);
  auto reference_msa = build_MSA_from_file(env->reference_file);
  auto tree = Tree(env->tree_file, reference_msa, env->model, query_msa);

  // tests
  auto pquerys = tree.place();
  vector<string> out;

  for (auto const &p : pquerys)
    out.push_back(pquery_to_jplace_string(p));

  EXPECT_EQ(out.size(), 2);

  // teardown

}
