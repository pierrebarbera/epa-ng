#pragma once

#include <type_traits>
#include <vector>

#include <cereal/types/vector.hpp>

#include "sample/Placement.hpp"
#include "seq/Sequence.hpp"

template< class Placement_Type >
class PQuery {

  public:
  using seqid_type     = size_t;
  using value_type     = Placement_Type;
  using iterator       = typename std::vector< value_type >::iterator;
  using const_iterator = typename std::vector< value_type >::const_iterator;

  PQuery() = default;

  template< typename T = Placement_Type,
            typename   = std::enable_if_t<
                std::is_same< T, Placement >::value > >
  PQuery( PQuery< Preplacement > const& other )
      : sequence_id_( other.sequence_id() )
      , placements_( other.size() )
  {
    auto const size = other.size();
    for( size_t i = 0; i < size; ++i ) {
      placements_[ i ] = Placement( other.at( i ) );
    }
  }
  PQuery( seqid_type const seq_id,
          std::string const& header )
      : sequence_id_( seq_id )
      , header_( header )
  {
  }
  PQuery( seqid_type const seq_id )
      : sequence_id_( seq_id )
  {
  }
  ~PQuery() = default;

  // move and copy semantics
  PQuery( PQuery const& other ) = default;
  PQuery( PQuery&& other )      = default;

  PQuery& operator=( PQuery const& other ) = default;
  PQuery& operator=( PQuery&& other ) = default;

  // needs to be in the header
  template< typename... Args >
  void emplace_back( Args&&... args )
  {
    placements_.emplace_back( std::forward< Args >( args )... );
  }

  // member access
  value_type& back() { return placements_.back(); }
  inline seqid_type sequence_id() const { return sequence_id_; }
  inline void sequence_id( seqid_type const seq_id ) { sequence_id_ = seq_id; }
  std::string const& header() const { return header_; }
  size_t size() const { return placements_.size(); }

  // manipulators
  void erase( iterator begin, iterator end ) { placements_.erase( begin, end ); }
  void resize( size_t size ) { return placements_.resize( size ); }
  inline void insert( iterator this_first, const_iterator begin, const_iterator end )
  {
    placements_.insert( this_first, begin, end );
  }

  inline void append( const_iterator begin, const_iterator end )
  {
    placements_.insert( placements_.end(), begin, end );
  }

  inline std::vector< value_type >& data()
  {
    return placements_;
  }

  // Iterator Compatibility
  iterator begin() { return placements_.begin(); }
  iterator end() { return placements_.end(); }
  const_iterator begin() const { return placements_.cbegin(); }
  const_iterator end() const { return placements_.cend(); }
  const_iterator cbegin() { return placements_.cbegin(); }
  const_iterator cend() { return placements_.cend(); }

  // Operator overloads
  value_type& operator[]( size_t const index ) { return placements_[ index ]; }
  value_type const& operator[]( size_t const index ) const { return placements_[ index ]; }
  value_type const& at( size_t const index ) const { return placements_.at( index ); }
  bool operator==( PQuery const&& other ) { return sequence_id_ == other.sequence_id_; }
  bool operator==( PQuery const& other ) { return sequence_id_ == other.sequence_id_; }

  // serialization
  template< class Archive >
  void serialize( Archive& ar ) { ar( sequence_id_, header_, placements_ ); }

  private:
  seqid_type sequence_id_ = 0;
  std::string header_;
  std::vector< value_type > placements_;
};
