#include "Epatest.hpp"
#include "check_equal.hpp"

#include <vector>

#include "core/pll/pllhead.hpp"
#include "core/raxml/Model.hpp"
#include "io/Binary.hpp"
#include "io/file_io.hpp"
#include "io/msa_reader.hpp"
#include "tree/Tree.hpp"
#include "util/Options.hpp"

using namespace std;

static void write_( Options const options )
{
  SKIP_CONFIG( options.memsave );
  // setup
  auto msa = build_MSA_from_file( env->reference_file, MSA_Info( env->reference_file ), options.premasking );
  raxml::Model model;

  Tree tree( env->tree_file, msa, model, options );

  // test
  dump_to_binary( tree, env->binary_file );
}

TEST( Binary, write )
{
  all_combinations( write_ );
}

static double loglh( pll_partition* partition, pll_unode_t* node )
{
  std::vector< unsigned int > param_indices( partition->rate_cats, 0 );
  return pll_compute_edge_loglikelihood( partition,
                                         node->clv_index,
                                         node->scaler_index,
                                         node->back->clv_index,
                                         node->back->scaler_index,
                                         node->pmatrix_index,
                                         &param_indices[ 0 ],
                                         nullptr );
}

static void read_( Options options )
{
  SKIP_CONFIG( options.memsave );
  // setup
  auto msa = build_MSA_from_file( env->reference_file, MSA_Info( env->reference_file ), options.premasking );
  raxml::Model model;
  Tree original_tree( env->tree_file, msa, model, options );
  dump_to_binary( original_tree, env->binary_file );

  // test
  Tree read_tree( env->binary_file, model, options );

  auto part      = original_tree.partition();
  auto read_part = read_tree.partition();

  // compare numbered jplace strings
  string original_nns( get_numbered_newick_string( original_tree.tree(), rtree_mapper() ) );
  string read_nns( get_numbered_newick_string( read_tree.tree(), rtree_mapper() ) );

  EXPECT_STREQ( original_nns.c_str(), read_nns.c_str() );
  // compare tree traversals
  ASSERT_EQ( original_tree.nums().nodes, read_tree.nums().nodes );
  vector< pll_unode_t* > original_nodes( original_tree.nums().nodes );
  vector< pll_unode_t* > read_nodes( read_tree.nums().nodes );
  unsigned int original_traversed, read_traversed;
  pll_utree_traverse( get_root( original_tree.tree() ),
                      PLL_TREE_TRAVERSE_POSTORDER,
                      cb_full_traversal,
                      &original_nodes[ 0 ],
                      &original_traversed );
  pll_utree_traverse( get_root( read_tree.tree() ),
                      PLL_TREE_TRAVERSE_POSTORDER,
                      cb_full_traversal,
                      &read_nodes[ 0 ],
                      &read_traversed );

  ASSERT_EQ( original_traversed, read_traversed );
  ASSERT_EQ( original_traversed, original_tree.nums().nodes );

  for( size_t i = 0; i < read_traversed; i++ ) {
    auto o = original_nodes[ i ];
    auto r = read_nodes[ i ];
    // printf("orig: %d back: %d\n", o->clv_index, o->back->clv_index);
    // printf("read: %d back: %d\n", r->clv_index, r->back->clv_index);
    check_equal( o, r );
    
    // ensure the read_tree has all CLVs, tipchars and scalers loaded
    read_tree.ensure_clv_loaded( r );
    if( r->next ) {
      read_tree.ensure_clv_loaded( r->next );
      read_tree.ensure_clv_loaded( r->next->next );
    }
  }

  // compare the partitions
  check_equal( *part, *read_part );

  // finally, compare some likelihoods
  for( size_t i = 0; i < read_traversed; i++ ) {
    auto o             = original_nodes[ i ];
    auto r             = read_nodes[ i ];
    auto original_logl = loglh( part, o );
    auto read_logl     = loglh( read_part, r );
    // printf("%f vs %f\n", original_logl, read_logl);
    EXPECT_DOUBLE_EQ( original_logl, read_logl );
  }
}

TEST( Binary, read )
{
  all_combinations( read_ );
}
