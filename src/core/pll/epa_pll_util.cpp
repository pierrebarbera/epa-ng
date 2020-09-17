#include "core/pll/epa_pll_util.hpp"

#include <algorithm>
#include <unordered_map>

#include "core/pll/optimize.hpp"
#include "core/pll/pll_util.hpp"
#include "set_manipulators.hpp"

void link_tree_msa( pll_utree_t* tree,
                    pll_partition_t* partition,
                    raxml::Model& model,
                    MSA const& msa,
                    unsigned int const num_tip_nodes )
{
  // assert( num_tip_nodes == msa.size() );

  // associate the sequences from the MSA file with the correct tips
  /* create a hash table of size num_tip_nodes */
  std::unordered_map< std::string, unsigned int > map; // mapping labels to tip clv indices

  /* populate the hash table with tree tip labels */
  for( size_t i = 0; i < num_tip_nodes; ++i ) {
    map[ tree->nodes[ i ]->label ] = tree->nodes[ i ]->clv_index;
  }

  /* find sequences in hash table and link them with the corresponding taxa */
  for( auto const& s : msa ) {
    auto map_value = map.find( s.header() );

    // failure tolerance: the MSA may also contain query sequences
    if( map_value == map.end() ) {
      continue;
      // throw std::runtime_error{std::string("Sequence with header does not appear in the tree: ") + s.header()};
    }

    auto clv_index = map_value->second;
    // associates the sequence with the tip by calculating the tips clv buffers
    pll_set_tip_states( partition, clv_index, model.charmap(), s.sequence().c_str() );

    // clear this ref taxon from the map
    map.erase( map_value );
  }

  // if map still contains entries, then some taxa were missing from the ref MSA. Fail and report!
  if( not map.empty() ) {
    LOG_ERR << "The reference Tree contained taxa that could not be found in the reference MSA: ";
    for( auto const& n : map ) {
      LOG_ERR << "Failed to find: " << n.first;
    }
    throw std::invalid_argument{ "Bad tree/ref msa combination (see errors above)" };
  }

  if( model.empirical_base_freqs() ) {
    compute_and_set_empirical_frequencies( partition, model );
  }

  raxml::assign( partition, model );
}

void precompute_clvs( pll_utree_t const* const tree,
                      pll_partition_t* partition,
                      Tree_Numbers const& nums )
{
  /* various buffers for creating a postorder traversal and operations structures */
  std::vector< unsigned int > param_indices( partition->rate_cats, 0 );
  std::vector< pll_unode_t* > travbuffer( nums.nodes );
  std::vector< double > branch_lengths( nums.branches );
  std::vector< unsigned int > matrix_indices( nums.branches );
  std::vector< pll_operation_t > operations( nums.nodes );

  auto const root = get_root( tree );

  utree_free_node_data( root );

  for( size_t i = 0; i < tree->tip_count; ++i ) {
    auto const node = tree->nodes[ i ];
    /* perform a partial postorder traversal of the unrooted tree  starting at the current tip
      and returning every node whose clv in the direction of the tip hasn't been calculated yet*/
    unsigned int traversal_size = 0;
    unsigned int num_matrices   = 0;
    unsigned int num_ops        = 0;
    if( pll_utree_traverse( node->back,
                            PLL_TREE_TRAVERSE_POSTORDER,
                            cb_partial_traversal,
                            &travbuffer[ 0 ],
                            &traversal_size )
        != PLL_SUCCESS ) {
      throw std::runtime_error{ "Function pll_unode_traverse() requires inner nodes as parameters" };
    }

    /* given the computed traversal descriptor, generate the operations
       structure, and the corresponding probability matrix indices that
       may need recomputing */
    pll_utree_create_operations( &travbuffer[ 0 ],
                                 traversal_size,
                                 &branch_lengths[ 0 ],
                                 &matrix_indices[ 0 ],
                                 &operations[ 0 ],
                                 &num_matrices,
                                 &num_ops );

    pll_update_prob_matrices( partition,
                              &param_indices[ 0 ], // use model 0
                              &matrix_indices[ 0 ], // matrices to update
                              &branch_lengths[ 0 ],
                              num_matrices ); // how many should be updated

    /* use the operations array to compute all num_ops inner CLVs. Operations
       will be carried out sequentially starting from operation 0 towrds num_ops-1 */
    pll_update_partials( partition, &operations[ 0 ], num_ops );
  }
  utree_free_node_data( root );
}

/* a callback function for traversing only down nodes with non-slotted CLVs
 * Assumes the node data pointer points to the partition!
 */
static int cb_traverse_unslotted( pll_unode_t* node )
{
  auto partition = static_cast< pll_partition_t* >( node->data );
  assert(partition);

  pll_clv_manager_t * clv_man = partition->clv_man;
  assert(clv_man);

  if( clv_man->slot_of_clvid[ node->clv_index ] != PLL_CLV_CLV_UNSLOTTED ) {
    // the CLV is slotted, so we:
    // 1) pin the clv in its slot
    clv_man->is_pinned[ node->clv_index ] = true;
    // 2) do not traverse down the subtree of this node. The way is shut.
    return 0;
  } else {
    // not slotted = traverse (and thus recompute)!
    return 1;
  }
}

/**
 * Compute all partials toward the branch defined by the given node, ending up
 * with valid CLVs at node->clv_index and node->back->clv_index.
 *
 * When using memory saver, saves recomputation by pinning already slotted CLVs
 * in place, such that they don't have to be computed. This is the intended use
 * of this function.
 *
 * @param tree          the overarching tree structure
 * @param nums          numbers associated with the tree
 * @param subtree_sizes the tree-associated subtree sizes information
 * @param node          node defining the target branch
 * @param partition     the partition
 */
void partial_compute_clvs(pll_utree_t* const tree,
                          Tree_Numbers const& nums,
                          unsigned int const * const subtree_sizes,
                          pll_unode_t* const node,
                          pll_partition_t* partition)
{
  // set vroot to selected node (and reverse this at the end)
  pll_unode_t* old_root = tree->vroot;
  tree->vroot = node;

  /*
   * Basic idea: set the data pointer of each node to point to the partition.
   * Then we do a largest-subtree-first traversal, with
   * a custom traversal callback that stops the traversal if a node's CLV is
   * currently slotted. The callback then also pins that CLV in place.
   */
  
  // go through all inner nodes, set data pointer to point to the partition
  size_t const nodes_count = tree->tip_count + tree->inner_count;
  for (size_t i = 0; i < nodes_count; ++i)
  {
    tree->nodes[ i ]->data = partition;
  }

  // traverse!
  /* various buffers for creating a postorder traversal and operations structures */
  std::vector< unsigned int > param_indices( partition->rate_cats, 0 );
  std::vector< pll_unode_t* > travbuffer( nums.nodes );
  std::vector< double > branch_lengths( nums.branches );
  std::vector< unsigned int > matrix_indices( nums.branches );
  std::vector< pll_operation_t > operations( nums.nodes );

  unsigned int traversal_size = 0;
  unsigned int num_matrices   = 0;
  unsigned int num_ops        = 0;

  if( pll_utree_traverse_lsf( tree,
                              subtree_sizes,
                              PLL_TREE_TRAVERSE_POSTORDER,
                              cb_traverse_unslotted,
                              &travbuffer[ 0 ],
                              &traversal_size )
      != PLL_SUCCESS ) {
    throw std::runtime_error{ std::string("pll_utree_traverse_lsf failed: ")
                              + pll_errmsg };
  }

  /* given the computed traversal descriptor, generate the operations
   structure, and the corresponding probability matrix indices that
   may need recomputing */
  pll_utree_create_operations( &travbuffer[ 0 ],
                               traversal_size,
                               &branch_lengths[ 0 ],
                               &matrix_indices[ 0 ],
                               &operations[ 0 ],
                               &num_matrices,
                               &num_ops );

  pll_update_prob_matrices( partition,
                            &param_indices[ 0 ], // use model 0
                            &matrix_indices[ 0 ], // matrices to update
                            &branch_lengths[ 0 ],
                            num_matrices ); // how many should be updated

  /* use the operations array to compute all num_ops inner CLVs. Operations
     will be carried out sequentially starting from operation 0 towrds num_ops-1 */
  pll_update_partials( partition, &operations[ 0 ], num_ops );



  // cleanup
  // unset the data pointers juuust incase some free() gets called
  for (size_t i = 0; i < nodes_count; ++i)
  {
    tree->nodes[ i ]->data = nullptr;
  }

  tree->vroot = old_root;
}

void split_combined_msa( MSA& source,
                         MSA& target,
                         Tree& tree )
{
  std::vector< pll_unode_t* > tip_nodes( tree.nums().tip_nodes );
  tip_nodes.assign( tree.tree()->nodes,
                    tree.tree()->nodes + tree.tree()->tip_count );

  auto falsegroup_begin = partition( source.begin(), source.end(),
                                     [tip_nodes = std::move( tip_nodes )]( Sequence const& em ) {
                                       return find( tip_nodes.begin(), tip_nodes.end(), em ) != tip_nodes.end();
                                     } );
  target.num_sites( source.num_sites() );
  target.move_sequences( falsegroup_begin, source.end() );
  source.erase( falsegroup_begin, source.end() );
}

bool operator==( pll_unode_t const* node, Sequence const& s )
{
  return s.header().compare( node->label ) == 0;
}

bool operator==( Sequence const& s, pll_unode_t const* node )
{
  return operator==( node, s );
}

raxml::Model get_model( pll_partition_t* partition )
{
  using namespace raxml;

  DataType seqtype = DataType::autodetect;

  if( partition->states == 4 ) {
    seqtype = DataType::dna;
  } else if( partition->states == 20 ) {
    seqtype = DataType::protein;
  } else {
    throw std::runtime_error{ "Couldn't determine sequence type from partition" };
  }

  Model model( seqtype );

  assign( model, partition );

  return model;
}
