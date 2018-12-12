#pragma once

#include <vector>
#include <utility>
#include <cassert>
#include <stdexcept>
#include <string>

class rtree_mapper
{
public:
  using map_type = std::vector<unsigned int>;

  rtree_mapper() = default;
  rtree_mapper( unsigned int root_edge,
                unsigned int proximal_edge,
                unsigned int distal_edge,
                double proximal_length,
                double distal_length,
                bool left )
    : utree_root_edge_(root_edge)
    , rtree_proximal_edge_(proximal_edge)
    , rtree_distal_edge_(distal_edge)
    , proximal_edge_length_(proximal_length)
    , distal_edge_length_(distal_length)
    , left_(left)
    {}

  ~rtree_mapper() = default;

  operator bool() const { return not map_.empty(); }

  map_type& map() { return map_; }
  map_type const& map() const { return map_; }
  void map( map_type& map ) { map_ = map; }
  void map( map_type&& map ) { map_ = std::move( map ); }

  bool uroot_is_left() const { return left_; }

  /*
    Transforms a branch_id and distal length of a placement in the utree
    into its equivalent in the rtree
   */
  auto in_rtree(  unsigned int branch_id,
                  double distal_length ) const
  {
    assert( not map_.empty() );
    assert( branch_id < map_.size() );

    if ( branch_id == utree_root_edge_ ) {
      // determine if it should land on proximal or distal rtree edge
      if ( distal_length > distal_edge_length_ ) {
        branch_id = rtree_proximal_edge_;
        // in this case we also need to adjust the distal length, as it is now on an edge
        // that had its direction to root flipped
        // how much length carries over?
        auto carryover = distal_length - distal_edge_length_;
        // thats how much we need to subtract from the proximal length
        distal_length = proximal_edge_length_ - carryover;
      } else {
        branch_id = rtree_distal_edge_;
      }

    } else {
      branch_id = map_[ branch_id ];
    }

    return std::make_pair( branch_id, distal_length );
  }

  unsigned int map_at( size_t const i ) const
  {
    if ( is_utree_root_edge( i ) ) {
      throw std::invalid_argument{std::string("Edge ") + std::to_string(i) + " is the root edge! Please handle separately"};
    }
    return map_.at(i);
  }

  auto distal_of_utree_root() const
  {
    return std::make_pair( rtree_distal_edge_, distal_edge_length_ );
  }

  auto proximal_of_utree_root() const
  {
    return std::make_pair( rtree_proximal_edge_, proximal_edge_length_ );
  }

  bool is_utree_root_edge( size_t const i ) const { return i == utree_root_edge_; }

private:
  /*
    needed to identify and deal with the special case of the root edge
   */
  unsigned int utree_root_edge_     = 0;
  unsigned int rtree_proximal_edge_ = 0;
  unsigned int rtree_distal_edge_   = 0;

  /*
    needed to correctly translate the distal lengths and to which rtree edge the placement goes
   */
  double proximal_edge_length_  = -1.0;
  double distal_edge_length_    = -1.0;

  bool left_ = true;

  /*
    maps from utree edge nums to rtree edge nums
   */
   map_type map_;
};