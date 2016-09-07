#include "Epatest.hpp"

#include <vector>

#include "src/Tree.hpp"
#include "src/Binary.hpp"
#include "src/Options.hpp"
#include "src/Model.hpp"
#include "src/file_io.hpp"

using namespace std;

void check_equal(pll_utree_t* a, pll_utree_t* b)
{
  EXPECT_DOUBLE_EQ(a->length, b->length);
  printf("AA %f is %f\n", a->length, b->length);

  if (!(a->next) || !(b->next))
  {
    ASSERT_FALSE(!(a->next) != !(b->next));
    return;
  }

  check_equal(a->back, b->back);

  a = a->next;
  b = b->next;
  EXPECT_DOUBLE_EQ(a->length, b->length);
  printf("BB %f is %f\n", a->length, b->length);

  check_equal(a->back, b->back);

  a = a->next;
  b = b->next;
  EXPECT_DOUBLE_EQ(a->length, b->length);
  printf("CC %f is %f\n", a->length, b->length);

  check_equal(a->back, b->back);

}


TEST(Binary, write)
{
  // setup
  auto msa = build_MSA_from_file(env->reference_file);
  Model model;
  Options options;
  Tree tree(env->tree_file, msa, model, options);

  // test
  dump_to_binary(tree, env->binary_file);

}

TEST(Binary, read)
{
  // setup
  auto msa = build_MSA_from_file(env->reference_file);
  Model model;
  double freqs[4] = {0.1,0.2,0.3,0.4};
  double alpha = 42.42;
  double subs[6] = {0.1,0.2,0.3,0.4,0.5,0.6};
  int symm[6] = {1,2,1,2,1,2};
  model.base_frequencies(freqs, 4);
  model.substitution_rates(subs, 6);
  model.symmetries(symm, 6);
  model.alpha(alpha);
  Options options;
  Tree original_tree(env->tree_file, msa, model, options);
  dump_to_binary(original_tree, env->binary_file);

  // test
  Tree read_tree(env->binary_file, options);

  // compare numbered jplace strings
  string original_nns(get_numbered_newick_string(original_tree.tree()));
  string read_nns(get_numbered_newick_string(read_tree.tree()->back));

  EXPECT_STREQ(original_nns.c_str(), read_nns.c_str());

  auto read_freqs = read_tree.model().base_frequencies();
  for (size_t i = 0; i < 4; i++)
    EXPECT_DOUBLE_EQ(freqs[i], read_freqs[i]);

  // auto read_alpha = read_tree.model().alpha();
  // EXPECT_DOUBLE_EQ(alpha, read_alpha);

  auto read_subs = read_tree.model().substitution_rates();
  for(size_t i = 0; i < 6; i++)
  {
    EXPECT_DOUBLE_EQ(subs[i], read_subs[i]);
  }

  // auto read_symm = read_tree.model().symmetries();
  // for(size_t i = 0; i < 6; i++)
  // {
  //   EXPECT_DOUBLE_EQ(symm[i], read_symm[i]);
  // }

  // auto utree = original_tree.tree();
  // auto read_utree = read_tree.tree();
  // if (!utree->next)
  //   utree = utree->back;
  // if (!read_utree->next)
  //   read_utree = read_utree->back;
  //
  // check_equal(utree, read_utree);

  // compare tree traversals
  ASSERT_EQ(original_tree.nums().nodes, read_tree.nums().nodes);
  vector<pll_utree_t *> original_nodes(original_tree.nums().nodes);
  vector<pll_utree_t *> read_nodes(read_tree.nums().nodes);
  unsigned int original_traversed, read_traversed;
  pll_utree_traverse(original_tree.tree(), cb_full_traversal, &original_nodes[0], &original_traversed);
  pll_utree_traverse(read_tree.tree(), cb_full_traversal, &read_nodes[0], &read_traversed);

  ASSERT_EQ(original_traversed, read_traversed);
  ASSERT_EQ(original_traversed, original_tree.nums().nodes);

  for (size_t i = 0; i < read_traversed; i++)
  {
    auto o = original_nodes[i];
    auto r = read_nodes[i];
    // printf("orig: %d back: %d\n", o->clv_index, o->back->clv_index);
    // printf("read: %d back: %d\n", r->clv_index, r->back->clv_index);
    EXPECT_EQ(o->clv_index, r->clv_index);
    EXPECT_DOUBLE_EQ(o->length, r->length);

    EXPECT_EQ(o->back->clv_index, r->back->clv_index);
    EXPECT_DOUBLE_EQ(o->back->length, r->back->length);

  }


  auto part = original_tree.partition();
  auto read_part = read_tree.partition();

  ASSERT_EQ(part->sites, read_part->sites);
  ASSERT_EQ(part->states, read_part->states);
  ASSERT_EQ(part->states_padded, read_part->states_padded);
  ASSERT_EQ(part->rate_cats, read_part->rate_cats);
  ASSERT_EQ(part->tips, read_part->tips);
  ASSERT_EQ(part->clv_buffers, read_part->clv_buffers);
  ASSERT_EQ(part->attributes, read_part->attributes);

  // compare tips
  if (read_part->attributes & PLL_ATTRIB_PATTERN_TIP)
    for (size_t i = 0; i < part->tips; i++)
    {
      pll_utree_t node;
      node.clv_index = i;
      node.scaler_index = 0;
      char* read_tipchars = (char*) read_tree.get_clv(&node);
      for (size_t j = 0; j < part->sites; j++)
      {
        EXPECT_EQ(part->tipchars[i][j], read_tipchars[j]);
      }
    }

  // compare clvs
  size_t start = (read_part->attributes & PLL_ATTRIB_PATTERN_TIP) ? part->tips : 0;

  for (size_t i = start; i < part->tips + part->clv_buffers; i++)
  {
    pll_utree_t node;
    node.clv_index = i;
    node.scaler_index = 0;
    double* read_clv = (double*) read_tree.get_clv(&node);
    for (size_t j = 0; j < part->sites + part->states_padded + part->rate_cats; j++)
    {
      ASSERT_DOUBLE_EQ(part->clv[i][j], read_clv[j]);
    }
  }

  // check scalers
  for (size_t i = 0; i < part->scale_buffers; i++)
  {
    pll_utree_t node;
    node.clv_index = 0;
    node.scaler_index = i;
    read_tree.get_clv(&node);
    for (size_t j = 0; j < part->sites; j++)
    {
      EXPECT_EQ(part->scale_buffer[i][j], read_part->scale_buffer[i][j]);
    }
  }
}
