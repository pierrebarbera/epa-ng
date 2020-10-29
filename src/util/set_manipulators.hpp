#pragma once

#include <limits>
#include <unordered_map>
#include <vector>

#include "core/Work.hpp"
#include "seq/MSA.hpp"
#include "util/Timer.hpp"
#include "util/logging.hpp"

/* Merging */
inline void merge( Timer<>& dest, Timer<> const& src )
{
  dest.insert( dest.end(), src.begin(), src.end() );
}

/**
 * Splits a <src> into <num_parts> number of equally sized <parts>.
 * <parts> is expected to be empty.
 *
 * REQUIRES T TO BE:
 *   Iterable
 *   has size()
 *   has insert()
 *
 * @param src    Container to split
 * @param parts     resulting parts vector
 * @param num_parts number of parts
 */
template< class T >
void split( T const& src,
            std::vector< T >& parts,
            unsigned int const num_parts )
{
  parts.clear();
  unsigned int chunk_size
      = ceil( src.size() / static_cast< double >( num_parts ) );
  auto move_begin = src.begin();

  while( move_begin < src.end() ) {
    auto move_end = std::min( move_begin + chunk_size, src.end() );
    parts.emplace_back();
    parts.back().insert( move_begin, move_end );
    move_begin = move_end;
  }
}

template< class T > void merge( T& dest, std::vector< T >&& parts )
{
  for( auto& p : parts ) {
    merge( dest, std::move( p ) );
  }
}

template< class T > void merge( T& dest, std::vector< T >& parts )
{
  for( auto const& p : parts ) {
    merge( dest, p );
  }
}

template< class T >
void merge( std::vector< T >& dest, std::vector< T > const& parts )
{
  dest.insert( std::end( dest ), std::begin( parts ), std::end( parts ) );
}
