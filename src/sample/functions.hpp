#pragma once

#include "sample/PQuery.hpp"
#include "sample/Placement.hpp"
#include "sample/Sample.hpp"
#include "util/logging.hpp"
#include "util/set_manipulators.hpp"

#include <limits>

template< class T > void sort_by_lwr( PQuery< T >& pq )
{
  sort( pq.begin(), pq.end(), []( T const& p_a, T const& p_b ) -> bool {
    return p_a.lwr() > p_b.lwr();
  } );
}

template< class T > void sort_by_logl( PQuery< T >& pq )
{
  std::sort( pq.begin(), pq.end(), []( T const& lhs, T const& rhs ) -> bool {
    return lhs.likelihood() > rhs.likelihood();
  } );
}

template< class T > void compute_and_set_lwr( Sample< T >& sample )
{
#ifdef __OMP
#pragma omp parallel for schedule( dynamic )
#endif
  for( size_t j = 0; j < sample.size(); ++j ) {
    auto& pq     = sample[ j ];
    double total = 0.0;

    // find the maximum
    auto max
        = std::max_element( pq.begin(),
                            pq.end(),
                            []( T const& lhs, T const& rhs ) {
                              return ( lhs.likelihood() < rhs.likelihood() );
                            } )
              ->likelihood();

    std::vector< double > exp_lh( pq.size() );
    // get the distances to the max
    for( size_t i = 0; i < pq.size(); ++i ) {
      exp_lh[ i ] = std::exp( pq[ i ].likelihood() - max );
      total += exp_lh[ i ];
    }

    // normalize the distances
    for( size_t i = 0; i < pq.size(); ++i ) {
      pq[ i ].lwr( exp_lh[ i ] / total );
    }
  }
}

template< class T >
typename PQuery< T >::iterator until_top_percent( PQuery< T >& pq, double const x )
{
  sort_by_lwr( pq );
  auto num_keep
      = static_cast< size_t >( ceil( x * static_cast< double >( pq.size() ) ) );
  auto iter = pq.begin();
  advance( iter, num_keep );
  return iter;
}

template< class T >
typename PQuery< T >::iterator until_accumulated_reached( PQuery< T >& pq,
                                                          double const thresh,
                                                          size_t const min,
                                                          size_t const max )
{
  sort_by_lwr( pq );

  double sum        = 0.0;
  size_t num_summed = 0;

  auto pq_iter = pq.begin();
  for( pq_iter = pq.begin(); num_summed < max and sum < thresh;
       ++pq_iter, ++num_summed ) {
    // sum up until threshold is passed. if we abort before it is passed, we
    // would have the possibility of empty lists
    sum += pq_iter->lwr();
  }

  auto to_add = distance( pq_iter, pq.begin() + min - 1 );
  if( to_add > 0 ) {
    advance( pq_iter, to_add );
  }

  return pq_iter;
}

template< class T >
typename PQuery< T >::iterator until_accumulated_reached( PQuery< T >& pq,
                                                          double const thresh )
{
  return until_accumulated_reached(
      pq, thresh, 1, std::numeric_limits< size_t >::max() );
}

template< class T >
void discard_bottom_x_percent( Sample< T >& sample, double const x )
{
  if( x < 0.0 || x > 1.0 ) {
    throw std::range_error { "x is not a percentage (outside of [0,1])" };
  }

#ifdef __OMP
#pragma omp parallel for schedule( dynamic )
#endif
  for( size_t i = 0; i < sample.size(); ++i ) {
    auto& pq        = sample[ i ];
    auto erase_iter = until_top_percent( pq, 1.0 - x );
    pq.erase( erase_iter, pq.end() );
  }
}

template< class T >
void discard_by_support_threshold( Sample< T >& sample,
                                   double const thresh,
                                   size_t const min = 1,
                                   size_t const max
                                   = std::numeric_limits< size_t >::max() )
{
  if( thresh < 0.0 or thresh > 1.0 ) {
    throw std::range_error {
      "thresh is not a valid likelihood weight ratio (outside of [0,1])"
    };
  }
  if( min < 1 ) {
    throw std::range_error { "Filter min cannot be smaller than 1!" };
  }

#ifdef __OMP
#pragma omp parallel for schedule( dynamic )
#endif
  for( size_t i = 0; i < sample.size(); ++i ) {
    auto& pq = sample[ i ];
    sort_by_lwr( pq );
    auto erase_iter
        = partition( pq.begin(), pq.end(), [ thresh ]( T& p ) -> bool {
            return ( p.lwr() > thresh );
          } );

    auto const num_kept
        = static_cast< size_t >( distance( pq.begin(), erase_iter ) );

    if( num_kept < min ) {
      auto const to_add = min - num_kept;
      std::advance( erase_iter, to_add );
    }

    if( max and num_kept > max ) {
      auto const to_remove = num_kept - max;
      std::advance( erase_iter, -to_remove );
    }

    pq.erase( erase_iter, pq.end() );
  }
}

template< class T >
void discard_by_accumulated_threshold( Sample< T >& sample,
                                       double const thresh,
                                       size_t const min = 1,
                                       size_t const max
                                       = std::numeric_limits< size_t >::max() )
{
  if( thresh < 0.0 || thresh > 1.0 ) {
    throw std::range_error {
      "thresh is not a valid likelihood weight ratio (outside of [0,1])"
    };
  }

  if( min < 1 ) {
    throw std::range_error { "Filter min cannot be smaller than 1!" };
  }

  if( min > max ) {
    throw std::range_error { "Filter min cannot be smaller than max!" };
  }

// accumulation and erasure phase
#ifdef __OMP
#pragma omp parallel for schedule( dynamic )
#endif
  for( size_t i = 0; i < sample.size(); ++i ) {
    auto& pq = sample[ i ];

    auto pq_iter = until_accumulated_reached( pq, thresh, min, max );

    pq.erase( pq_iter, pq.end() );
  }
}

template< class T > void filter( Sample< T >& sample, Options const& options )
{
  if( options.acc_threshold ) {
    LOG_DBG << "Filtering output by accumulated threshold: "
            << options.support_threshold << std::endl;
    discard_by_accumulated_threshold( sample,
                                      options.support_threshold,
                                      options.filter_min,
                                      options.filter_max );
  } else {
    LOG_DBG << "Filtering output placements below threshold: "
            << options.support_threshold << std::endl;
    discard_by_support_threshold( sample,
                                  options.support_threshold,
                                  options.filter_min,
                                  options.filter_max );
  }
}

/**
 * collapses PQuerys with the same ID inside a Sample into one
 */
template< class T > void collapse( Sample< T >& sample )
{
  auto const invalid
      = std::numeric_limits< typename PQuery< T >::seqid_type >::max();

  std::unordered_map< size_t, std::vector< size_t > > collapse_set;

  // build map of all pqueries
  for( size_t i = 0; i < sample.size(); ++i ) {
    auto const& pq = sample[ i ];
    collapse_set[ pq.sequence_id() ].emplace_back( i );
  }

  // find all cases of puplicate entries and merge hem
  for( auto& pair : collapse_set ) {
    auto pqlist = pair.second;
    // duplicate!
    if( pqlist.size() > 1 ) {
      // move entries from duplicate to original
      auto& dest = sample[ pqlist[ 0 ] ].data();
      for( size_t i = 1; i < pqlist.size(); ++i ) {
        auto& src = sample[ pqlist[ i ] ].data();
        dest.reserve( dest.size() + src.size() );
        std::move(
            std::begin( src ), std::end( src ), std::back_inserter( dest ) );
        // mark invalid in original sample
        sample[ pqlist[ i ] ].sequence_id( invalid );
      }
    }
  }

  // clear the original sample of invalid pqueries
  sample.erase( std::remove_if( std::begin( sample ),
                                std::end( sample ),
                                [ invalid = invalid ]( auto& e ) {
                                  return e.sequence_id() == invalid;
                                } ),
                std::end( sample ) );
}

/**
 * special split function that Splits samples in buckets according to the global sequence ID
 * of their PQueries. The goal is to have them split such that each aggregate node gets their
 * correct set of sequence results (even if that part is empty, which constitutes a null-message)
 *
 * @param
 * @param
 * @param
 */
template< class T >
void split( Sample< T > const& src,
            std::vector< Sample< T > >& parts,
            unsigned int const num_parts )
{
  parts.clear();
  // ensure that there are actually as many parts as specified. We want empty parts to enable null messages
  parts.resize( num_parts );

  for( auto& pq : src ) {
    auto const bucket = pq.sequence_id() % num_parts;
    parts[ bucket ].push_back( pq );
  }
}

/**
  Merges a Sample <src> into a Sample <dest>. Leaves <src> intact.
*/
template< class T >
void merge( Sample< T >& dest, Sample< T > const& src )
{
  // merge in every source pquery...
  for( auto const& pquery : src ) {
    // ... by checking if its sequence already exists in destination
    auto input_iter = find( dest.begin(), dest.end(), pquery );
    // if not, create a record
    if( input_iter == dest.end() ) {
      dest.emplace_back( pquery.sequence_id(), pquery.header() );
      input_iter = --( dest.end() );
    }
    // then concat their vectors
    input_iter->insert( input_iter->end(), pquery.begin(), pquery.end() );
  }
}

/**
  Merges a Sample <src> into a Sample <dest>. src here is an rvalue,
  and thus the elements are moved instead of copied
*/
template< class T >
void merge( Sample< T >& dest, Sample< T >&& src )
{
  for( auto& pquery : src ) {
    // create new record
    dest.emplace_back( std::move( pquery ) );
  }
}
