#pragma once

#include "core/pll/pllhead.hpp"
#include "io/file_io.hpp"
#include "seq/MSA.hpp"
#include "seq/MSA_Info.hpp"
#include "tree/Tiny_Tree.hpp"
#include "tree/Tree.hpp"
#include "util/Options.hpp"

inline void check_equal_subnode( pll_unode_t const* const a,
                                 pll_unode_t const* const b )
{
  EXPECT_STREQ( a->label, a->label );
  EXPECT_DOUBLE_EQ( a->length, b->length );
  EXPECT_EQ( a->node_index, b->node_index );
  EXPECT_EQ( a->clv_index, b->clv_index );
  EXPECT_EQ( a->scaler_index, b->scaler_index );
  EXPECT_EQ( a->pmatrix_index, b->pmatrix_index );
}

inline void check_equal( pll_unode_t const* const a,
                         pll_unode_t const* const b )
{
  check_equal_subnode( a, b );

  if( a->next ) {
    ASSERT_TRUE( b->next != nullptr );
    check_equal_subnode( a->next, b->next );
    check_equal_subnode( a->next->next, b->next->next );
  }
}

inline void check_equal( pll_utree_t const& lhs, pll_utree_t const& rhs )
{
  ASSERT_EQ( lhs.tip_count, rhs.tip_count );
  ASSERT_EQ( lhs.inner_count, rhs.inner_count );
  ASSERT_EQ( lhs.edge_count, rhs.edge_count );
  EXPECT_EQ( lhs.binary, rhs.binary );

  for( size_t i = 0; i < lhs.tip_count + lhs.inner_count; ++i ) {
    check_equal( lhs.nodes[ i ], rhs.nodes[ i ] );
  }
}

// inline std::vector< unsigned int > scaler_to_clv_map( Tree const& tree )
// {
//   auto const num_scalers = tree.partition()->scale_buffers;
//   std::vector< unsigned int > map( num_scalers );

//   std::vector< pll_unode_t* > travbuffer( tree.nums().nodes );
//   unsigned int trav_size = 0;
//   pll_utree_traverse( get_root( tree.tree() ),
//                       PLL_TREE_TRAVERSE_POSTORDER,
//                       full_trav,
//                       &travbuffer[ 0 ],
//                       &trav_size );

//   for( auto& n : travbuffer ) {
//     if( n->scaler_index != PLL_SCALE_BUFFER_NONE ) {
//       map[ n->scaler_index ] = n->clv_index;
//     }
//   }

//   return map;
// }

inline void check_equal( pll_partition_t const& lhs,
                         pll_partition_t const& rhs)
{
  // TODO this function ignores asc_bias
  // TODO this function ignores p_inv

  EXPECT_EQ( lhs.attributes, rhs.attributes );
  auto const attributes = lhs.attributes;
  ASSERT_FALSE( attributes & PLL_ATTRIB_LIMIT_MEMORY );

  ASSERT_EQ( lhs.clv_buffers, rhs.clv_buffers );
  ASSERT_EQ( lhs.tips, rhs.tips );
  ASSERT_EQ( lhs.nodes, rhs.nodes );
  ASSERT_EQ( lhs.pattern_weight_sum, rhs.pattern_weight_sum );
  ASSERT_EQ( lhs.rate_matrices, rhs.rate_matrices );
  ASSERT_EQ( lhs.prob_matrices, rhs.prob_matrices );
  ASSERT_EQ( lhs.sites, rhs.sites );
  ASSERT_EQ( lhs.rate_cats, rhs.rate_cats );
  ASSERT_EQ( lhs.states, rhs.states );
  ASSERT_EQ( lhs.scale_buffers, rhs.scale_buffers );

  auto const clv_buffers   = lhs.clv_buffers;
  auto const scale_buffers = lhs.scale_buffers;
  auto const tips          = lhs.tips;
  auto const sites         = lhs.sites;
  auto const states        = lhs.states;

  // compare CLVs
  size_t const start
      = ( lhs.attributes & PLL_ATTRIB_PATTERN_TIP ) ? lhs.tips : 0;
  for( size_t clv_id = start; clv_id < clv_buffers + tips; ++clv_id ) {
    double const* lhs_clv = lhs.clv[ clv_id ];
    double const* rhs_clv = rhs.clv[ clv_id ];
    auto const clv_size   = pll_get_clv_size( &lhs, clv_id );
    ASSERT_EQ( clv_size, pll_get_clv_size( &rhs, clv_id ) );
    for( size_t i = 0; i < clv_size; ++i ) {
      EXPECT_DOUBLE_EQ( lhs_clv[ i ], rhs_clv[ i ] );
    }
  }

  for( size_t i = 0; i < states; i++ ) {
    EXPECT_DOUBLE_EQ( lhs.frequencies[ 0 ][ i ], rhs.frequencies[ 0 ][ i ] );
  }

  for( size_t i = 0; i < 6; i++ ) {
    EXPECT_DOUBLE_EQ( lhs.subst_params[ 0 ][ i ], rhs.subst_params[ 0 ][ i ] );
  }

  for( size_t i = 0; i < lhs.rate_cats; i++ ) {
    EXPECT_DOUBLE_EQ( lhs.rates[ i ], rhs.rates[ i ] );
  }

  for( size_t i = 0; i < lhs.rate_cats; i++ ) {
    EXPECT_DOUBLE_EQ( lhs.rate_weights[ i ], rhs.rate_weights[ i ] );
  }

  // compare tips
  if( attributes & PLL_ATTRIB_PATTERN_TIP ) {
    if( lhs.tipchars and rhs.tipchars ) {
      for( size_t tip_id = 0; tip_id < lhs.tips; tip_id++ ) {
        for( size_t site = 0; site < sites; site++ ) {
          if( lhs.tipchars[ tip_id ] and rhs.tipchars[ tip_id ] ) {
            auto const lhs_chr = lhs.tipchars[ tip_id ][ site ];
            auto const rhs_chr = rhs.tipchars[ tip_id ][ site ];

            ASSERT_EQ( lhs_chr, rhs_chr )
                << "(tip_id=" << tip_id << ", site=" << site << ")";
          } else {
            ASSERT_EQ( lhs.tipchars[ tip_id ], rhs.tipchars[ tip_id ] );
          }
        }
      }
    } else {
      ASSERT_EQ( lhs.tipchars, rhs.tipchars );
    }
  }

  // check scalers
  for( size_t i = 0; i < scale_buffers; i++ ) {

    auto const scaler_size = pll_get_sites_number( &lhs, i );
    ASSERT_EQ( scaler_size, pll_get_sites_number( &rhs, i ) );
    if( lhs.scale_buffer[ i ] and rhs.scale_buffer[ i ] ) {
      for( size_t j = 0; j < scaler_size; j++ ) {
        // printf("%u v %u\n",lhs.scale_buffer[i][j], rhs.scale_buffer[i][j] );
        EXPECT_EQ( lhs.scale_buffer[ i ][ j ], rhs.scale_buffer[ i ][ j ] );
      }
    } else {
      ASSERT_EQ( lhs.scale_buffer[ i ], rhs.scale_buffer[ i ] );
    }
  }

  // check repeats specific structures
  auto const lhs_rep = lhs.repeats;
  auto const rhs_rep = rhs.repeats;
  if( lhs_rep ) {
    ASSERT_TRUE( rhs_rep );
    for( size_t scaler_id = 0; scaler_id < scale_buffers; ++scaler_id ) {
      EXPECT_EQ( lhs_rep->perscale_ids[ scaler_id ],
                 rhs_rep->perscale_ids[ scaler_id ] );
    }
    for( size_t clv_id = 0; clv_id < clv_buffers + tips; ++clv_id ) {
      EXPECT_EQ( lhs_rep->pernode_ids[ clv_id ],
                 rhs_rep->pernode_ids[ clv_id ] );
      EXPECT_EQ( lhs_rep->pernode_allocated_clvs[ clv_id ],
                 rhs_rep->pernode_allocated_clvs[ clv_id ] );

      for( size_t j = 0; j < sites; ++j ) {
        EXPECT_EQ( lhs_rep->pernode_site_id[ clv_id ][ j ],
                   rhs_rep->pernode_site_id[ clv_id ][ j ] );
      }

      for( size_t j = 0; j < pll_get_sites_number( &lhs, clv_id ); ++j ) {
        EXPECT_EQ( lhs_rep->pernode_id_site[ clv_id ][ j ],
                   rhs_rep->pernode_id_site[ clv_id ][ j ] );
      }
    }
  }
}

inline void check_equal( Tiny_Tree const& lhs, Tiny_Tree const& rhs )
{
  EXPECT_EQ( lhs.branch_id(), rhs.branch_id() );
  EXPECT_DOUBLE_EQ( lhs.distal_length(), rhs.distal_length() );
  EXPECT_DOUBLE_EQ( lhs.proximal_length(), rhs.proximal_length() );
  EXPECT_DOUBLE_EQ( lhs.pendant_length(), rhs.pendant_length() );

  check_equal( *lhs.utree(), *rhs.utree() );
  check_equal( *lhs.partition(), *rhs.partition() );
}
