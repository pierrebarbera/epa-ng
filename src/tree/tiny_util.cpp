#include "tree/tiny_util.hpp"

#include <cstring>
#include <type_traits>

#include "core/pll/pll_util.hpp"
#include "core/pll/error.hpp"

constexpr unsigned int proximal_clv_index        = 4;
constexpr unsigned int inner_clv_index           = 3;
constexpr unsigned int new_tip_clv_index         = 1;
constexpr unsigned int distal_clv_index_if_tip   = 2;
constexpr unsigned int distal_clv_index_if_inner = 5;

template< class T >
static void alloc_and_copy( T*& dest,
                            T const* const src,
                            size_t const size,
                            size_t const alignment = 0ul )
{
  if( dest != nullptr ) {
    free( dest );
  }

  if( alignment ) {
    dest = static_cast< T* >(
        pll_aligned_alloc( size * sizeof( T ), alignment ) );
  } else {
    dest = static_cast< T* >( calloc( size, sizeof( T ) ) );
  }

  if( dest == nullptr ) {
    throw std::runtime_error { "Can't alloc memory." };
  }

  assert( dest != nullptr );

  std::memcpy( dest,
               src,
               size * sizeof( T ) );
}

static void shallow_copy_clv(pll_partition_t* dest_part,
                           pll_unode_t const* const dest_node,
                           pll_partition_t const* const src_part,
                           pll_unode_t const* const src_node)
{
  pll_aligned_free( dest_part->clv[ dest_node->clv_index ] );
  assert( src_part->clv[ src_node->clv_index ] );
  dest_part->clv[ dest_node->clv_index ]
      = src_part->clv[ src_node->clv_index ];
}


static void deep_copy_clv( pll_partition_t* dest_part,
                           pll_unode_t const* const dest_node,
                           pll_partition_t const* const src_part,
                           pll_unode_t const* const src_node )
{
  /**
   * TODO:
   * - this function should probably cover the repeats case as well
   */

  // clv size in doubles (gets correct size wrt repeats also)
  auto const clv_size = pll_get_clv_size( src_part,
                                          src_node->clv_index );

  // get the correct clv pointer (depends on memory mode)
  double const* src_clv = pll_get_clv_reading( src_part,
                                               src_node->clv_index );

  if( src_part == nullptr ) {
    throw std::runtime_error { "Failed to fetch source clv." };
  }

  double* dest_clv = dest_part->clv[ dest_node->clv_index ];

  // assert( dest_clv != nullptr );

  if( dest_clv != nullptr ) {
    // copy the clv buffer over
    std::memcpy( dest_clv,
                 src_clv,
                 clv_size * sizeof( double ) );
  } else {
    alloc_and_copy( dest_part->clv[ dest_node->clv_index ],
                    src_clv,
                    clv_size,
                    src_part->alignment );
  }
}

static void deep_copy_scaler( pll_partition_t* dest_part,
                              pll_unode_t* dest_node,
                              pll_partition_t const* const src_part,
                              pll_unode_t const* const src_node )
{
  if( src_node->scaler_index != PLL_SCALE_BUFFER_NONE
      and src_part->scale_buffer[ src_node->scaler_index ] != nullptr ) {

    // sites in the associated CLV (correct number if repeats)
    auto const sites = pll_get_sites_number( src_part, src_node->clv_index );

    auto const sites_alloc = src_part->asc_additional_sites + sites;
    auto const scaler_size = ( src_part->attributes & PLL_ATTRIB_RATE_SCALERS )
        ? sites_alloc * src_part->rate_cats
        : sites_alloc;

    alloc_and_copy( dest_part->scale_buffer[ dest_node->scaler_index ],
                    src_part->scale_buffer[ src_node->scaler_index ],
                    scaler_size );
  }
}

static void deep_copy_repeats( pll_partition_t* dest_part,
                               pll_unode_t* dest_node,
                               pll_partition_t const* const src_part,
                               pll_unode_t const* const src_node )
{
  // copy size info
  if( src_node->scaler_index != PLL_SCALE_BUFFER_NONE ) {
    dest_part->repeats->perscale_ids[ dest_node->scaler_index ]
        = src_part->repeats->perscale_ids[ src_node->scaler_index ];
  }

  dest_part->repeats->pernode_ids[ dest_node->clv_index ]
      = src_part->repeats->pernode_ids[ src_node->clv_index ];
  dest_part->repeats->pernode_allocated_clvs[ dest_node->clv_index ]
      = src_part->repeats->pernode_allocated_clvs[ src_node->clv_index ];

  if( src_part->repeats->pernode_ids[ src_node->clv_index ] ) {
    auto const size = pll_get_sites_number( const_cast< pll_partition_t* >( src_part ),
                                            src_node->clv_index );

    alloc_and_copy( dest_part->repeats->pernode_site_id[ dest_node->clv_index ],
                    src_part->repeats->pernode_site_id[ src_node->clv_index ],
                    src_part->sites );
    alloc_and_copy( dest_part->repeats->pernode_id_site[ dest_node->clv_index ],
                    src_part->repeats->pernode_id_site[ src_node->clv_index ],
                    size );
  }
}

pll_partition_t* make_tiny_partition( pll_partition_t const* const old_partition,
                                      pll_utree_t const* tree,
                                      pll_unode_t const* const old_proximal,
                                      pll_unode_t const* const old_distal,
                                      bool const tip_tip_case,
                                      bool const deep_copy_clvs )
{
  /**
    As we work with PLL_PATTERN_TIP functionality, special care has to be taken
    in regards to the node and partition structure: PLL assumes that any node
    with clv index < number of tips is in fact a real tip, that is a tip that
    uses a character array instead of a real clv. Here we need to set up the
    node/partition to fool pll: the tips that actually contain CLVs copied over
    from the reference node have their index set to greater than number of tips.
    This results in a acceptable amount of wasted memory that is never used
    (num_sites * bytes * number of clv-tips)
  */
  assert( old_partition );
  assert( tree );
  assert( old_proximal );
  assert( old_distal );

  // ensure: never shallow copy from memory managed partitions.
  assert( not( ( old_partition->attributes & PLL_ATTRIB_LIMIT_MEMORY )
               and not deep_copy_clvs ) );

  unsigned int attributes = old_partition->attributes;

  // unset memory saver mode for the tiny partition
  attributes &= ~PLL_ATTRIB_LIMIT_MEMORY;

  bool const use_tipchars = attributes & PLL_ATTRIB_PATTERN_TIP;

  // tip_inner case: both reference nodes are inner nodes
  // tip tip case: one for the "proximal" clv tip
  unsigned int const num_clv_tips = tip_tip_case ? 1 : 2;

  auto proximal = tree->nodes[ 0 ];
  auto distal   = tree->nodes[ 1 ];

  pll_partition_t* tiny = pll_partition_create(
      3, // tips
      1 + num_clv_tips, // extra clvs
      old_partition->states,
      old_partition->sites,
      old_partition->rate_matrices,
      3, // number of prob. matrices (one per possible unique branch length)
      old_partition->rate_cats,
      3, // number of scale buffers (one per possible inner node)
      attributes );

  assert( tiny );

  unsigned int i;
  free( tiny->rates );
  tiny->rates = old_partition->rates;
  if( tiny->subst_params ) {
    for( i = 0; i < tiny->rate_matrices; ++i ) {
      pll_aligned_free( tiny->subst_params[ i ] );
    }
  }
  free( tiny->subst_params );
  tiny->subst_params = old_partition->subst_params;
  if( tiny->frequencies ) {
    for( i = 0; i < tiny->rate_matrices; ++i ) {
      pll_aligned_free( tiny->frequencies[ i ] );
    }
  }
  free( tiny->frequencies );
  tiny->frequencies = old_partition->frequencies;
  if( tiny->eigenvecs ) {
    for( i = 0; i < tiny->rate_matrices; ++i ) {
      pll_aligned_free( tiny->eigenvecs[ i ] );
    }
  }
  free( tiny->eigenvecs );
  tiny->eigenvecs = old_partition->eigenvecs;
  if( tiny->inv_eigenvecs ) {
    for( i = 0; i < tiny->rate_matrices; ++i ) {
      pll_aligned_free( tiny->inv_eigenvecs[ i ] );
    }
  }
  free( tiny->inv_eigenvecs );
  tiny->inv_eigenvecs = old_partition->inv_eigenvecs;
  if( tiny->eigenvals ) {
    for( i = 0; i < tiny->rate_matrices; ++i ) {
      pll_aligned_free( tiny->eigenvals[ i ] );
    }
  }
  free( tiny->eigenvals );
  tiny->eigenvals = old_partition->eigenvals;

  if( tiny->prop_invar ) {
    free( tiny->prop_invar );
  }
  tiny->prop_invar = old_partition->prop_invar;

  if( tiny->invariant ) {
    free( tiny->invariant );
  }
  tiny->invariant = old_partition->invariant;

  free( tiny->eigen_decomp_valid );
  tiny->eigen_decomp_valid = old_partition->eigen_decomp_valid;
  if( tiny->pattern_weights ) {
    free( tiny->pattern_weights );
  }
  tiny->pattern_weights = old_partition->pattern_weights;

  // copy the repeats structures
  if( old_partition->repeats ) {
    // then do the per-clv stuff, but only for the two relevant clv
    deep_copy_repeats( tiny,
                       proximal,
                       old_partition,
                       old_proximal );

    deep_copy_repeats( tiny,
                       distal,
                       old_partition,
                       old_distal );

    pll_resize_repeats_lookup( tiny, tiny->sites * tiny->states );
  }

  if( use_tipchars ) {
    // alloc the tipchars array
    assert( tiny->tipchars == nullptr );
    tiny->tipchars = static_cast< unsigned char** >(
        calloc( tiny->tips, sizeof( unsigned char* ) ) );

    // hardcopy the charmap buffers
    alloc_and_copy( tiny->charmap, old_partition->charmap, PLL_ASCII_SIZE );
    alloc_and_copy( tiny->tipmap, old_partition->tipmap, PLL_ASCII_SIZE );
    alloc_and_copy( tiny->ttlookup,
                    old_partition->ttlookup,
                    size_of_ttlookup( old_partition ),
                    old_partition->alignment );
    tiny->maxstates = old_partition->maxstates;

    // alloc the memory for the pendant/query tip
    auto const sites_alloc
        = tiny->asc_bias_alloc ? tiny->sites + tiny->states : tiny->sites;
    tiny->tipchars[ new_tip_clv_index ] = static_cast< unsigned char* >(
        calloc( sites_alloc, sizeof( unsigned char ) ) );
  }

  // handle the distal
  if( tip_tip_case and use_tipchars ) {
    tiny->tipchars[ distal->clv_index ]
        = old_partition->tipchars[ old_distal->clv_index ];
  } else {
    if( deep_copy_clvs ) {
      deep_copy_clv( tiny, distal, old_partition, old_distal );
    } else {
      shallow_copy_clv( tiny, distal, old_partition, old_distal );
    }
  }

  // handle the proximal
  if( deep_copy_clvs ) {
    deep_copy_clv( tiny, proximal, old_partition, old_proximal );
  } else {
    shallow_copy_clv( tiny, proximal, old_partition, old_proximal );
  }

  // deep copy scalers
  deep_copy_scaler( tiny,
                    proximal,
                    old_partition,
                    old_proximal );

  deep_copy_scaler( tiny,
                    distal,
                    old_partition,
                    old_distal );


  return tiny;
}

void tiny_partition_destroy( pll_partition_t* partition,
                             bool const deep_copy_clvs )
{
  if( partition ) {
    // unset shallow copied things
    // TODO recheck thath these can be reused
    partition->rates              = nullptr;
    partition->subst_params       = nullptr;
    partition->frequencies        = nullptr;
    partition->eigenvecs          = nullptr;
    partition->inv_eigenvecs      = nullptr;
    partition->eigenvals          = nullptr;
    partition->prop_invar         = nullptr;
    partition->invariant          = nullptr;
    partition->eigen_decomp_valid = nullptr;
    partition->pattern_weights    = nullptr;

    if( not deep_copy_clvs ) {
      partition->clv[ proximal_clv_index ] = nullptr;
    }

    bool const distal_is_tip    = partition->clv_buffers == 3 ? false : true;
    bool const pattern_tip_mode = partition->attributes & PLL_ATTRIB_PATTERN_TIP;

    if( distal_is_tip ) {
      if( pattern_tip_mode ) {
        partition->tipchars[ distal_clv_index_if_tip ] = nullptr;
      } else if ( not deep_copy_clvs ) {
          partition->clv[ distal_clv_index_if_tip ] = nullptr;
      }
    } else if( not deep_copy_clvs ) {
      partition->clv[ distal_clv_index_if_inner ] = nullptr;
    }

    pll_partition_destroy( partition );
  }
}

/*
  destructor that assumes partition was shallow copied
 */
void tiny_partition_destroy_shallow( pll_partition_t* partition )
{
  tiny_partition_destroy( partition, false );
}

/*
  destructor that assumes partition was deep copied
 */
void tiny_partition_destroy_deep( pll_partition_t* partition )
{
  tiny_partition_destroy( partition, true );
}

pll_utree_t* make_tiny_tree_structure( pll_unode_t const* old_proximal,
                                       pll_unode_t const* old_distal,
                                       double const original_brlen,
                                       bool const tip_tip_case )
{
  unsigned int const inner_scaler_index    = 1;
  unsigned int const proximal_scaler_index = 0;
  unsigned int const distal_scaler_index   = 2;

  /**
    As we work with PLL_PATTERN_TIP functionality, special care has to be taken
    in regards to the tree and partition structure: PLL assumes that any node
    with clv index < number of tips is in fact a real tip, that is a tip that
    uses a character array instead of a real clv. Here we need to set up the
    tree/partition to fool pll: the tips that actually contain CLVs copied over
    from the reference tree have their index set to greater than number of tips.
    This results in a acceptable amount of wasted memory that is never used
    (num_sites * bytes
    * number of clv-tips)
  */
  // if tip-inner case
  unsigned int distal_clv_index
      = ( tip_tip_case ) ? distal_clv_index_if_tip : distal_clv_index_if_inner;

  auto inner
      = static_cast< pll_unode_t* >( calloc( 1, sizeof( pll_unode_t ) ) );
  inner->next
      = static_cast< pll_unode_t* >( calloc( 1, sizeof( pll_unode_t ) ) );
  inner->next->next
      = static_cast< pll_unode_t* >( calloc( 1, sizeof( pll_unode_t ) ) );
  inner->next->next->next = inner;

  auto new_tip
      = static_cast< pll_unode_t* >( calloc( 1, sizeof( pll_unode_t ) ) );
  auto proximal
      = static_cast< pll_unode_t* >( calloc( 1, sizeof( pll_unode_t ) ) );
  auto distal
      = static_cast< pll_unode_t* >( calloc( 1, sizeof( pll_unode_t ) ) );

  // connect the nodes to each other
  inner->back             = new_tip;
  new_tip->back           = inner;
  inner->next->back       = distal;
  distal->back            = inner->next;
  inner->next->next->back = proximal;
  proximal->back          = inner->next->next;

  // set up clv indices
  inner->clv_index             = inner_clv_index;
  inner->next->clv_index       = inner_clv_index;
  inner->next->next->clv_index = inner_clv_index;
  proximal->clv_index          = proximal_clv_index;
  distal->clv_index            = distal_clv_index;
  new_tip->clv_index           = new_tip_clv_index;

  // set up scaler indices
  new_tip->scaler_index           = PLL_SCALE_BUFFER_NONE;
  inner->scaler_index             = inner_scaler_index;
  inner->next->scaler_index       = inner_scaler_index;
  inner->next->next->scaler_index = inner_scaler_index;
  proximal->scaler_index
      = ( old_proximal->scaler_index == PLL_SCALE_BUFFER_NONE )
      ? PLL_SCALE_BUFFER_NONE
      : proximal_scaler_index;
  distal->scaler_index = ( old_distal->scaler_index == PLL_SCALE_BUFFER_NONE )
      ? PLL_SCALE_BUFFER_NONE
      : distal_scaler_index;

  reset_triplet_lengths( inner, nullptr, original_brlen );

  auto tree = static_cast< pll_utree_t* >( calloc( 1, sizeof( pll_utree_t ) ) );

  tree->binary = true;

  tree->edge_count  = 3;
  tree->tip_count   = 3;
  tree->inner_count = 1;

  tree->nodes = static_cast< pll_unode_t** >(
      calloc( tree->inner_count + tree->tip_count,
              sizeof( pll_unode_t* ) ) );

  tree->nodes[ 0 ] = proximal;
  tree->nodes[ 1 ] = distal;
  tree->nodes[ 2 ] = new_tip;
  tree->nodes[ 3 ] = inner;
  tree->vroot      = inner;

  return tree;
}
