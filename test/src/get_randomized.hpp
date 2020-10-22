#pragma once

#include <random>
#include <cmath>

#include "core/Work.hpp"

constexpr int SEED = 42;

// these functions would be great to templetize like get_randomized< T >,
// but thats a yak for another shaving session

inline Work get_randomized_Work( size_t const num_branches,
                                 size_t const num_queries,
                                 double const draw_fraction )
{
  Work work;

  std::mt19937 gen;
  gen.seed( SEED );
  std::uniform_int_distribution<> branch_dist( 0, num_branches );
  std::uniform_int_distribution<> query_dist( 0, num_queries );

  auto const num_drawn_branches
      = static_cast< size_t >( std::ceil( num_branches * draw_fraction ) );
  auto const num_drawn_queries = static_cast< size_t >(
      std::ceil( ( num_queries * draw_fraction ) / num_drawn_branches ) );

  for( size_t i = 0; i < num_drawn_branches; ++i ) {
    for( size_t j = 0; j < num_drawn_queries; ++j ) {
      work.add( branch_dist( gen ), query_dist( gen ) );
    }
  }

  return work;
}
