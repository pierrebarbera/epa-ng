#pragma once

#include <string>
#include <vector>

class Sequence {
  public:
  Sequence()  = default;
  ~Sequence() = default;
  Sequence( std::string header, std::string sequence )
      : sequence_( sequence )
  {
    header_.push_back( header );
  }
  Sequence( Sequence const& s ) = default;
  Sequence( Sequence&& s )      = default;

  // operator overloads
  Sequence& operator=( Sequence const& s ) = default;
  Sequence& operator=( Sequence&& s ) = default;
  bool operator==( Sequence const& other ) { return sequence_.compare( other.sequence() ) == 0; }
  bool operator==( Sequence const& other ) const { return sequence_.compare( other.sequence() ) == 0; }
  // TODO doesn't merge in the full list (very tailored to the collapse func)
  void merge( Sequence const& other ) { header_.push_back( other.header() ); }

  // member access
  std::string const& header() const { return header_.front(); }
  std::vector< std::string > const& header_list() const { return header_; }
  std::string const& sequence() const { return sequence_; }

  private:
  std::vector< std::string > header_;
  std::string sequence_;
};
