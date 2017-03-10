#include "Epatest.hpp"

#include "src/pllhead.hpp"
#include "src/file_io.hpp"
#include "src/Tree_Numbers.hpp"
#include "src/Model.hpp"
#include "src/Tiny_Tree.hpp"
#include "src/Tree.hpp"
#include "src/Binary.hpp"
#include "src/MSA.hpp"
#include "src/Range.hpp"
#include "src/pll_util.hpp"
#include "src/epa_pll_util.hpp"
#include "src/Lookup_Store.hpp"

#include <tuple>
#include <limits>

using namespace std;

TEST(Tiny_Tree, place_heuristic)
{
  // buildup
  MSA msa = build_MSA_from_file(env->reference_file);
  MSA queries = build_MSA_from_file(env->query_file);
  Options options;

  auto ref_tree = Tree(env->tree_file, msa, env->model, options);

  shared_ptr<Lookup_Store> lu_ptr(new Lookup_Store(ref_tree.nums().branches, ref_tree.partition()->states));

  auto tree = ref_tree.tree();

  // tests
  Tiny_Tree tt(tree, 0, ref_tree, false, options, lu_ptr);

  for (auto const &x : queries)
  {
    auto place = tt.place(x);
    EXPECT_NE(place.likelihood(), 0.0);
    // printf("lk: %f\n", place.likelihood());
    ASSERT_NE (place.likelihood(), -numeric_limits<double>::infinity());
    EXPECT_NE(place.distal_length(), 0.0);
    EXPECT_NE(place.pendant_length(), 0.0);
  }
  // teardown
}

TEST(Tiny_Tree, place_BLO)
{
  // buildup
  MSA msa = build_MSA_from_file(env->reference_file);
  MSA queries = build_MSA_from_file(env->query_file);


  Options options;
  auto ref_tree = Tree(env->tree_file, msa, env->model, options);

  shared_ptr<Lookup_Store> lu_ptr(new Lookup_Store(ref_tree.nums().branches, ref_tree.partition()->states));

  auto tree = ref_tree.tree();

  // tests
  Tiny_Tree tt(tree, 0, ref_tree, true, options, lu_ptr);

  for (auto const &x : queries)
  {
    auto place = tt.place(x);
    EXPECT_NE(place.likelihood(), 0.0);
    EXPECT_NE(place.distal_length(), 0.0);
    EXPECT_NE(place.pendant_length(), 0.0);
  }
  // teardown
}

TEST(Tiny_Tree, place_from_binary)
{
  // setup
  auto msa = build_MSA_from_file(env->reference_file);
  auto queries = build_MSA_from_file(env->query_file);
  Model model;
  Options options;

  Tree original_tree(env->tree_file, msa, model, options);
  dump_to_binary(original_tree, env->binary_file);
  Tree read_tree(env->binary_file, model, options);
  string invocation("./this --is -a test");
  shared_ptr<Lookup_Store> lu_ptr(new Lookup_Store(original_tree.nums().branches, original_tree.partition()->states));

  ASSERT_EQ(original_tree.nums().branches, read_tree.nums().branches);

  vector<pll_utree_t *> original_branches(original_tree.nums().branches);
  vector<pll_utree_t *> read_branches(read_tree.nums().branches);
  auto original_traversed  = utree_query_branches(original_tree.tree(), &original_branches[0]);
  auto read_traversed  = utree_query_branches(read_tree.tree(), &read_branches[0]);

  ASSERT_EQ(original_traversed, read_traversed);
  ASSERT_EQ(original_traversed, original_tree.nums().branches);

  // test
  for (size_t i = 0; i < original_traversed; i++)
  {
    Tiny_Tree original_tiny(original_branches[i], 0, original_tree, false, options, lu_ptr);
    Tiny_Tree read_tiny(read_branches[i], 0, read_tree, false, options, lu_ptr);
    for(auto& seq : queries)
    {
      auto original_place = original_tiny.place(seq);
      auto read_place = read_tiny.place(seq);

      ASSERT_DOUBLE_EQ(original_place.likelihood(), read_place.likelihood());
      // printf("%f vs %f\n", original_place.likelihood(), read_place.likelihood());
    }
  }
  // teardown
}
