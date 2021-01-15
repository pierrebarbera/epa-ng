#include "core/pll/epa_pll_util.hpp"

#include <algorithm>
#include <unordered_map>
#include <vector>

#include "core/pll/optimize.hpp"
#include "core/pll/pll_util.hpp"
#include "util/Options.hpp"
#include "util/logging.hpp"

/**
 * makes a vector mapping from internal node_indexes to their associated branch
 * IDs (where the branch ids are consistent with the jplace standard)
 */
std::vector< unsigned int > get_branch_ids( pll_utree_t const* const tree )
{
  // size of the result vector will be number of tips + number of internal
  // nodes (which are triplets)
  std::vector< unsigned int > result( 3 * tree->inner_count + tree->tip_count );

  std::vector< pll_unode_t* > branches( tree->edge_count );
  auto traversed = utree_query_branches( tree,
                                         &branches[ 0 ] );

  for( unsigned int branch_id = 0; branch_id < traversed; ++branch_id ) {
    auto const& node = branches[ branch_id ];
    // both node and node back obviously are associated with the same branch
    result[ node->node_index ] = result[ node->back->node_index ] = branch_id;
  }

  return result;
}

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
    handle_pll_failure( not pll_utree_traverse( node->back,
                                                PLL_TREE_TRAVERSE_POSTORDER,
                                                cb_partial_traversal,
                                                &travbuffer[ 0 ],
                                                &traversal_size ),
                        "utree traverse failed during clv precomputation" );

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

using clv_cost_pair
    = std::pair< decltype( pll_unode_t::clv_index ), unsigned int >;

typedef struct {
  unsigned int const* subtree_sizes;
  pll_partition_t* partition;
  std::vector< clv_cost_pair > pin_list;
} partial_subtree_sizes_data;

/**
 * Returns true if the vector contains the element of the corresponding clv index.
 * If the vector does contain it, remove it from the vector.
 */
static bool contains_consume( std::vector< clv_cost_pair >& vec,
                              decltype( pll_unode_t::clv_index ) clv_index )
{
  auto iter = std::find_if(
      vec.begin(), vec.end(), [ clv_index ]( clv_cost_pair const& elem ) {
        return elem.first == clv_index;
      } );

  if( iter != vec.end() ) {
    vec.erase( iter );
    return true;
  } else {
    return false;
  }
}

/**
 * callback to identify candiates for pinning.
 *
 * Assumes the node data pointer points to the partition!
 *
 * @param  node current node
 * @return      whether to keep traversing
 */
static int cb_traverse_ident_pinnable( pll_unode_t* node )
{
  auto cb_data = static_cast< partial_subtree_sizes_data* >( node->data );
  assert( cb_data );

  auto partition = cb_data->partition;
  assert( partition );

  pll_clv_manager_t* clv_man = partition->clv_man;
  assert( clv_man );

  if( clv_man->slot_of_clvid[ node->clv_index ] != PLL_CLV_CLV_UNSLOTTED ) {
    // the CLV is slotted, add it to the list and stop traversing
    auto const cost = cb_data->subtree_sizes[ node->node_index ];
    cb_data->pin_list.emplace_back( node->clv_index, cost );

    return 0;
  } else {
    // not slotted = traverse
    return 1;
  }
}

/* a callback function for traversing only down nodes with non-slotted CLVs
 * Assumes the node data pointer points to the partition!
 */
static int cb_traverse_unslotted( pll_unode_t* node )
{
  auto cb_data = static_cast< partial_subtree_sizes_data* >( node->data );
  assert( cb_data );

  auto partition = cb_data->partition;
  assert( partition );

  pll_clv_manager_t* clv_man = partition->clv_man;
  assert( clv_man );

  if( clv_man->slot_of_clvid[ node->clv_index ] != PLL_CLV_CLV_UNSLOTTED
      and contains_consume( cb_data->pin_list, node->clv_index ) ) {
    // the CLV is slotted AND we have marked it to be pinned
    // 1) pin the clv in its slot

    auto retval = pll_pin_clv( partition, node->clv_index );
    assert( retval );
    (void)retval;

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
void partial_compute_clvs( pll_utree_t* const tree,
                           Tree_Numbers const& nums,
                           unsigned int const* const subtree_sizes,
                           pll_unode_t* const node,
                           pll_partition_t* partition )
{
  // set vroot to selected node (and reverse this at the end)
  pll_unode_t* old_root = tree->vroot;
  tree->vroot           = node->next ? node : node->back;

  /*
   * Basic idea: set the data pointer of each node to point to a custom struct
   * containing pointers to the partition, the subtree_sizes array, and a list
   * of clv candidates for pinning. Then we do two largest-subtree-first
   * traversals, one to fill the candidates list, and one to pin the selected
   * candidates (and abort the traversal at the pinned nodes/clvs). Inbetween
   * these two iterations, we make a selection from the list of those we could
   * pin, by keeping only as many as we can afford and pinning just the most
   * valuable ones.
   *
   * Addendum: as this function is intended to be called for ever changing
   * vroots (as we progress through the tree according to some traversal)
   * and we pin in place CLVs that are already computed, we will sometimes
   * end up in a state where we forced clvs to be slotted such that,
   * effectively, it looks like we "traversed" in a non-lsf order, and thus run
   * out of slots. The solution to this is to check every time we would pin an
   * existing CLV whether it would leave us with sufficient clv slots to finish
   * the computation, and not pin it if it doesn't
   */

  // go through all inner nodes, set data pointer to point to the partition
  partial_subtree_sizes_data cb_data;
  cb_data.subtree_sizes = subtree_sizes;
  cb_data.partition     = partition;

  size_t const nodes_count = tree->tip_count + tree->inner_count;
  for( size_t i = 0; i < nodes_count; ++i ) {
    auto cur_node = tree->nodes[ i ];
    cur_node->data = &cb_data;
    if( cur_node->next ) {
      cur_node->next->data = cur_node->next->next->data = &cb_data;
    }
  }

  // set up for the traversals
  std::vector< pll_unode_t* > pretrav_buf( nums.nodes );
  unsigned int pretrav_size = 0;

  // First, traverse once to identify slotted CLVs that could be pinned
  handle_pll_failure( not pll_utree_traverse_lsf( tree,
                                                  subtree_sizes,
                                                  PLL_TREE_TRAVERSE_POSTORDER,
                                                  cb_traverse_ident_pinnable,
                                                  &pretrav_buf[ 0 ],
                                                  &pretrav_size ),
                      "pll_utree_traverse_lsf failed." );

  // then, keep only the highest cost entries of the list of clvs we could pin

  // figure out how many extra slots we can pin
  size_t const lower_free_slot_bound = ceil( log2( tree->tip_count ) + 2 );
  auto const pinning_budget
      = partition->clv_man->slottable_size - lower_free_slot_bound;

  // sort ccandidates by cost
  std::sort( std::begin( cb_data.pin_list ),
             std::end( cb_data.pin_list ),
             [](clv_cost_pair& lhs, clv_cost_pair& rhs) {
                return lhs.second > rhs.second;
             } );

  // discard the lower ones
  auto candidate_iter = std::begin( cb_data.pin_list );
  std::advance( candidate_iter, pinning_budget );
  cb_data.pin_list.erase( candidate_iter,
                                    std::end( cb_data.pin_list ) );

  /* various buffers for creating a postorder traversal and operations structures */
  std::vector< unsigned int > param_indices( partition->rate_cats, 0 );
  std::vector< pll_unode_t* > travbuffer( nums.nodes );
  std::vector< double > branch_lengths( nums.branches );
  std::vector< unsigned int > matrix_indices( nums.branches );
  std::vector< pll_operation_t > operations( nums.nodes );

  unsigned int traversal_size = 0;
  unsigned int num_matrices   = 0;
  unsigned int num_ops        = 0;

  handle_pll_failure( not pll_utree_traverse_lsf( tree,
                                                  subtree_sizes,
                                                  PLL_TREE_TRAVERSE_POSTORDER,
                                                  cb_traverse_unslotted,
                                                  &travbuffer[ 0 ],
                                                  &traversal_size ),
                      "pll_utree_traverse_lsf failed." );

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

  // preparation for any next iteration: unpin the CLVs at node and node->back
  // (especially needed if next iteration has different root direction, as
  // otherwise one of these CLV will never get unpinned by update_partials)
  pll_unpin_clv( partition, node->clv_index );
  pll_unpin_clv( partition, node->back->clv_index );

  // assert that none of the CLVs are pinned
  assert( partition->clv_man->num_pinned == 0 );
  assert( [ & ]() {
    for( auto clv_index = partition->clv_man->addressable_begin;
         clv_index < partition->clv_man->addressable_end;
         ++clv_index ) {
      if( partition->clv_man->is_pinned[ clv_index ] ) {
        return false;
      }
    }
    return true;
  }() );

  // cleanup
  // unset the data pointers juuust incase some free() gets called
  for( size_t i = 0; i < nodes_count; ++i ) {
    auto cur_node = tree->nodes[ i ];
    cur_node->data = nullptr;
    if( cur_node->next ){
      cur_node->next->data = cur_node->next->next->data = nullptr;
    }
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

static unsigned int simd_autodetect()
{
  if( PLL_STAT( avx2_present ) )
    return PLL_ATTRIB_ARCH_AVX2;
  else if( PLL_STAT( avx_present ) )
    return PLL_ATTRIB_ARCH_AVX;
  else if( PLL_STAT( sse3_present ) )
    return PLL_ATTRIB_ARCH_SSE;
  else
    return PLL_ATTRIB_ARCH_CPU;
}

pll_partition_t* make_partition( raxml::Model const& model,
                                 pll_utree_t* const tree,
                                 Tree_Numbers const& nums,
                                 int const num_sites,
                                 Options const& options,
                                 Memory_Config const& memsave )
{
  assert( nums.tip_nodes ); // nums must have been initialized correctly

  auto attributes = simd_autodetect();

  if( ( options.scaling == Options::Numerical_Scaling::kOn )
      or ( ( options.scaling == Options::Numerical_Scaling::kAuto )
           and nums.large_tree() ) ) {
    attributes |= PLL_ATTRIB_RATE_SCALERS;
  }

  if( options.repeats ) {
    attributes |= PLL_ATTRIB_SITE_REPEATS;
  } else {
    attributes |= PLL_ATTRIB_PATTERN_TIP;
  }

  if( memsave ) {
    if( options.repeats ) {
      throw std::runtime_error{ "Repeats + memsave not supported" };
    }
    attributes |= PLL_ATTRIB_PATTERN_TIP;
    attributes |= PLL_ATTRIB_LIMIT_MEMORY;
  }

  auto partition = pll_partition_create(
      nums.tip_nodes,
      nums.inner_nodes
          * 3, // number of extra clv buffers: 3 for every direction on the node
      model.num_states(),
      num_sites,
      1,
      nums.branches,
      model.num_ratecats(),
      ( nums.inner_nodes * 3 ) + nums.tip_nodes, /* number of scaler buffers */
      attributes );

  handle_pll_failure(
    not partition, "Could not create partition (make_partition)");

  if( memsave ) {
    assert( memsave.structs.subtree_sizes() );

    handle_pll_failure( not pll_clv_manager_init(
                            partition, memsave.clv_slots, NULL, NULL, NULL ),
                        "Could not initialize CLV manager." );

    handle_pll_failure(
        not pll_clv_manager_MRC_strategy_init(
            partition->clv_man, tree, memsave.structs.subtree_sizes() ),
        "Could not initialize CLV manager replacement strategy" );
  }
  return partition;
}
