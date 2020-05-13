#pragma once

#include <utility>
#include <vector>

#include "seq/Sequence.hpp"

class MSA {
  public:
  using container_type = std::vector< Sequence >;
  using iterator       = container_type::iterator;
  using const_iterator = container_type::const_iterator;

  explicit MSA( size_t const num_sites )
      : num_sites_( num_sites ){};
  MSA()
      : num_sites_( 0 ){};
  ~MSA() = default;

  void move_sequences( iterator begin, iterator end );
  void append( std::string const& header, std::string const& sequence );
  void erase( iterator begin, iterator end ) { sequence_list_.erase( begin, end ); }
  void clear() { sequence_list_.clear(); }

  static void swap( MSA& a, MSA& b )
  {
    std::swap( a.num_sites_, b.num_sites_ );
    std::swap( a.sequence_list_, b.sequence_list_ );
  }

  // getters
  size_t size() const { return sequence_list_.size(); }
  size_t num_sites() const { return num_sites_; }
  Sequence const& operator[]( size_t const i ) const { return sequence_list_[ i ]; }

  // setters
  void num_sites( size_t const sites ) { num_sites_ = sites; }

  //Iterator Compatability
  iterator begin() { return sequence_list_.begin(); }
  iterator end() { return sequence_list_.end(); }
  const_iterator begin() const { return sequence_list_.cbegin(); }
  const_iterator end() const { return sequence_list_.cend(); }
  const_iterator cbegin() const { return sequence_list_.cbegin(); }
  const_iterator cend() const { return sequence_list_.cend(); }

  private:
  size_t num_sites_;
  container_type sequence_list_;
};

namespace std {
void swap( MSA& a, MSA& b );
}
