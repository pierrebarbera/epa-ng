#include "Epatest.hpp"

#include <vector>

#include "src/Tree.hpp"
#include "src/Binary.hpp"
#include "src/Options.hpp"
#include "src/Model.hpp"
#include "src/file_io.hpp"

using namespace std;

void check_equal(pll_unode_t* a, pll_unode_t* b)
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
  options.opt_branches = options.opt_model = true;

  Tree tree(env->tree_file, msa, model, options);

  // test
  dump_to_binary(tree, env->binary_file);

}

static int full_trav(pll_unode_t*)
{
  return 1;
}

static auto create_scaler_to_clv_map(Tree& tree)
{
  const auto num_scalers = tree.partition()->scale_buffers;
  std::vector<unsigned int> map(num_scalers);

  std::vector<pll_unode_t*> travbuffer(tree.nums().nodes);
  unsigned int trav_size = 0;
  pll_utree_traverse( get_root(tree.tree()),
                      PLL_TREE_TRAVERSE_POSTORDER,
                      full_trav, 
                      &travbuffer[0], 
                      &trav_size);

  for (auto& n : travbuffer) {
    if (n->scaler_index != PLL_SCALE_BUFFER_NONE) {
      map[n->scaler_index] = n->clv_index;
    }
  }

  return map;
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
  Tree read_tree(env->binary_file, model, options);

  // compare numbered jplace strings
  string original_nns(get_numbered_newick_string(original_tree.tree()));
  string read_nns(get_numbered_newick_string(read_tree.tree()));

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
  vector<pll_unode_t *> original_nodes(original_tree.nums().nodes);
  vector<pll_unode_t *> read_nodes(read_tree.nums().nodes);
  unsigned int original_traversed, read_traversed;
  pll_utree_traverse( get_root(original_tree.tree()),
                      PLL_TREE_TRAVERSE_POSTORDER,
                      cb_full_traversal, 
                      &original_nodes[0], 
                      &original_traversed);
  pll_utree_traverse( get_root(read_tree.tree()),
                      PLL_TREE_TRAVERSE_POSTORDER,
                      cb_full_traversal,
                      &read_nodes[0],
                      &read_traversed);

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
      pll_unode_t node;
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
    pll_unode_t node;
    node.clv_index = i;
    node.scaler_index = 0;
    double* read_clv = (double*) read_tree.get_clv(&node);
    const auto clv_size = pll_get_clv_size(part, i);
    for (size_t j = 0; j < clv_size; j++)
    {
      ASSERT_DOUBLE_EQ(part->clv[i][j], read_clv[j]);
    }
  }

  const auto stoc = create_scaler_to_clv_map(read_tree);

  // check scalers
  for (size_t i = 0; i < part->scale_buffers; i++)
  {
    pll_unode_t node;
    node.clv_index = 0;
    node.scaler_index = i;
    read_tree.get_clv(&node);
    const auto scaler_size =  pll_get_sites_number(read_part, stoc[i]);
    for (size_t j = 0; j < scaler_size; j++)
    {
      EXPECT_EQ(part->scale_buffer[i][j], read_part->scale_buffer[i][j]);
    }
  }

  // check repeats specific structures
  const auto rep = part->repeats; 
  const auto read_rep = read_part->repeats; 
  if (rep) {
    for (size_t i = 0; i < part->scale_buffers; ++i) {
      EXPECT_EQ(rep->perscale_max_id[i], read_rep->perscale_max_id[i]);
    }
    for (size_t i = 0; i < part->clv_buffers + part->tips; ++i) {
      EXPECT_EQ(rep->pernode_max_id[i], read_rep->pernode_max_id[i]);
      EXPECT_EQ(rep->pernode_allocated_clvs[i], read_rep->pernode_allocated_clvs[i]);

      for (size_t j = 0; j < part->sites; ++j) {
        EXPECT_EQ(rep->pernode_site_id[i][j], read_rep->pernode_site_id[i][j]);
      }

      for (size_t j = 0; j < pll_get_sites_number(part, i); ++j) {
        EXPECT_EQ(rep->pernode_id_site[i][j], read_rep->pernode_id_site[i][j]);
      }
    }
  }
}
