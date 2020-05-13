#pragma once

#include <algorithm>
#include <string>
#include <type_traits>
#include <vector>

#include <cereal/types/base_class.hpp>
#include <cereal/types/string.hpp>
#include <cereal/types/vector.hpp>

#include "pipeline/Token.hpp"
#include "sample/PQuery.hpp"

template< class Placement_Type = Placement >
class Sample : public Token {
  public:
  using value_type     = PQuery< Placement_Type >;
  using iterator       = typename std::vector< value_type >::iterator;
  using const_iterator = typename std::vector< value_type >::const_iterator;

  Sample() = default;

  template< typename T = Placement_Type,
            typename   = std::enable_if_t<
                std::is_same< T, Placement >::value > >
  Sample( Sample< Slim_Placement > const& other )
      : pquerys_( other.size() )
      , newick_( other.newick() )
  {
    auto const size = other.size();
    for( size_t i = 0; i < size; ++i ) {
      pquerys_[ i ] = PQuery< Placement >( other.at( i ) );
    }
  }
  explicit Sample( size_t const size )
  {
    pquerys_.reserve( size );
    for( size_t i = 0; i < size; ++i ) {
      pquerys_.emplace_back( i );
    }
  }

  Sample( size_t const size, size_t const depth )
      : Sample( size )
  {
    for( size_t i = 0; i < pquerys_.size(); ++i ) {
      pquerys_[ i ].resize( depth );
    }
  }

  explicit Sample( std::string const& newick )
      : newick_( newick )
  {
  }
  ~Sample() = default;

  // member access
  value_type& back() { return pquerys_.back(); }
  unsigned int size() const { return pquerys_.size(); }
  std::string const& newick() const { return newick_; }
  void clear() { pquerys_.clear(); }
  void push_back( value_type&& pq ) { pquerys_.push_back( pq ); }
  void push_back( value_type& pq ) { pquerys_.push_back( pq ); }
  void push_back( value_type const& pq ) { pquerys_.push_back( pq ); }
  void erase( iterator begin, iterator end ) { pquerys_.erase( begin, end ); }

  // needs to be in the header
  template< typename... Args >
  void emplace_back( Args&&... args ) { pquerys_.emplace_back( std::forward< Args >( args )... ); }
  template< class InputIt >
  void insert( InputIt first, InputIt last ) { pquerys_.insert( pquerys_.end(), first, last ); }

  template< typename... Args >
  void add_placement( size_t const seq_id,
                      std::string const& label,
                      Args&&... args )
  {
    // if seq_id in pquerys_
    auto iter = std::end( pquerys_ );
    if( ( iter = std::find( std::begin( pquerys_ ), std::end( pquerys_ ), value_type( seq_id ) ) )
        != std::end( pquerys_ ) ) {
      iter->emplace_back( std::forward< Args >( args )... );
    } else {
      pquerys_.emplace_back( seq_id, label );
      pquerys_.back().emplace_back( std::forward< Args >( args )... );
    }
  }

  size_t add_pquery( size_t const seq_id,
                     std::string const& label )
  {
    // Create a new pquery and return its id in the vector.
    pquerys_.emplace_back( seq_id, label );
    return pquerys_.size() - 1;
  }

  // Iterator Compatibility
  iterator begin() { return pquerys_.begin(); }
  iterator end() { return pquerys_.end(); }
  const_iterator begin() const { return pquerys_.cbegin(); }
  const_iterator end() const { return pquerys_.cend(); }
  const_iterator cbegin() { return pquerys_.cbegin(); }
  const_iterator cend() { return pquerys_.cend(); }

  // Operator overloads
  value_type& operator[]( size_t const index ) { return pquerys_[ index ]; }
  value_type const& at( size_t const index ) const { return pquerys_[ index ]; }

  // serialization
  template< class Archive >
  void serialize( Archive& ar )
  {
    ar( *static_cast< Token* >( this ), pquerys_, newick_ );
  }

  private:
  std::vector< value_type > pquerys_;
  std::string newick_;
};
