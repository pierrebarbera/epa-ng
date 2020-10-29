#pragma once

#include "core/Work.hpp"
#include "sample/Sample.hpp"
#include "sample/functions.hpp"
#include "util/Options.hpp"

#include <algorithm>
#include <type_traits>
#include <utility>

#ifdef __OMP
#include <omp.h>
#endif

static inline size_t get_num_threads( Options const& options )
{
#ifdef __OMP
  size_t const num_threads
      = options.num_threads ? options.num_threads : omp_get_max_threads();
  omp_set_num_threads( num_threads );
#else
  (void)options;
  size_t const num_threads = 1;
#endif
  return num_threads;
}

static inline size_t get_thread_id()
{
#ifdef __OMP
  return omp_get_thread_num();
#else
  return 0;
#endif
}

template< class T, typename F >
static Work
heuristic_( Sample< T >& sample, Options const& options, F& filterstop )
{
  Work result;
  compute_and_set_lwr( sample );

  auto const num_threads = get_num_threads( options );

  std::vector< Work > workvec( num_threads );

#ifdef __OMP
#pragma omp parallel for schedule( dynamic )
#endif
  for( size_t i = 0; i < sample.size(); ++i ) {
    auto& pq       = sample[ i ];
    auto const tid = get_thread_id();

    auto end = filterstop( pq, options.prescoring_threshold );

    for( auto iter = pq.begin(); iter != end; ++iter ) {
      workvec[ tid ].add( iter->branch_id(), pq.sequence_id() );
    }
  }
  merge( result, workvec );
  return result;
}

template< typename... Ts >
using until_type
    = decltype( until_accumulated_reached( std::declval< Ts >()... ) )( Ts... );

template< class T >
Work dynamic_heuristic( Sample< T >& sample, Options const& options )
{
  return heuristic_< T, until_type< PQuery< T >&, double const > >(
      sample, options, until_accumulated_reached );
}

template< class T >
Work fixed_heuristic( Sample< T >& sample, Options const& options )
{
  return heuristic_< T, decltype( until_top_percent< T > ) >(
      sample, options, until_top_percent );
}

template< class T >
Work baseball_heuristic( Sample< T >& sample, Options const& options )
{
  Work result;

  auto const num_threads = get_num_threads( options );

  // strike_box: logl delta, keep placements within this many logl units from
  // the best
  double const strike_box = 3;
  // max_strikes: number of additional branches to add after strike box is full
  size_t const max_strikes = 6;
  // max_pitches: absolute maximum of candidates to select
  size_t const max_pitches = 40;

  std::vector< Work > workvec( num_threads );
#ifdef __OMP
#pragma omp parallel for schedule( dynamic )
#endif
  for( size_t i = 0; i < sample.size(); ++i ) {
    auto& pq       = sample[ i ];
    auto const tid = get_thread_id();

    assert( pq.size() );
    // sort placements by likelihood (descending)
    sort_by_logl( pq );
    // keep any placements that are within strike box of the best
    double const best_logl = pq[ 0 ].likelihood();
    double const thresh    = best_logl - strike_box;
    // get first element not within strike box
    auto keep_iter
        = std::find_if( pq.begin(), pq.end(), [ thresh ]( auto const& p ) {
            return ( p.likelihood() < thresh );
          } );

    auto const hits = std::distance( pq.begin(), keep_iter );

    // ensure we keep no more than max_pitches
    size_t to_add = std::min( max_pitches - hits, max_strikes );

    std::advance( keep_iter, to_add );

    for( auto iter = pq.begin(); iter != keep_iter; ++iter ) {
      workvec[ tid ].add( iter->branch_id(), pq.sequence_id() );
    }
  }
  merge( result, workvec );
  return result;
}

template< class T >
Work apply_heuristic( Sample< T >& sample, Options const& options )
{
  if( options.baseball ) {
    return baseball_heuristic( sample, options );
  } else if( options.prescoring_by_percentage ) {
    return fixed_heuristic( sample, options );
  } else {
    return dynamic_heuristic( sample, options );
  }
}
