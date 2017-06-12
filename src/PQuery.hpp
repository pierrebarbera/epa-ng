#pragma once

#include <vector>
#include <cereal/types/vector.hpp>

#include "Sequence.hpp"
#include "Placement.hpp"

class PQuery {
public:
  using value_type      = Placement;
  using iterator        = std::vector<value_type>::iterator;
  using const_iterator  = std::vector<value_type>::const_iterator;

  PQuery() = default;
  PQuery (const unsigned int seq_id, const unsigned int size)
    : sequence_id_(seq_id), placements_(size) 
  { }
  PQuery (const unsigned int seq_id)
    : sequence_id_(seq_id) 
  { }
  ~PQuery() = default;

  // move and copy semantics
  PQuery(PQuery const& other) = default;
  PQuery(PQuery&& other)      = default;

  PQuery& operator= (PQuery const& other) = default;
  PQuery& operator= (PQuery && other)     = default;

  // needs to be in the header
  template<typename ...Args> void emplace_back(Args && ...args)
  {
    placements_.emplace_back(std::forward<Args>(args)...);
  }

  // member access
  Placement& back() { return placements_.back(); }
  unsigned int sequence_id() const { return sequence_id_; }
  double entropy() const { return entropy_; }
  void entropy(const double e) { entropy_ = e; }
  unsigned int size() const { return placements_.size(); }
  void erase(iterator begin, iterator end) { placements_.erase(begin, end); }

  void insert(iterator this_first, const_iterator begin, const_iterator end)
  { placements_.insert(this_first, begin, end); }

  // Iterator Compatibility
  iterator begin() { return placements_.begin(); }
  iterator end() { return placements_.end(); }
  const_iterator begin() const { return placements_.cbegin(); }
  const_iterator end() const { return placements_.cend(); }
  const_iterator cbegin() { return placements_.cbegin(); }
  const_iterator cend() { return placements_.cend(); }

  // Operator overloads
  Placement& operator[] (const unsigned int index) { return placements_[index]; }
  bool operator==(const PQuery&& other) { return sequence_id_ == other.sequence_id_; }
  bool operator==(const PQuery& other) { return sequence_id_ == other.sequence_id_; }

  // serialization
  template<class Archive>
  void serialize(Archive& ar) { ar( sequence_id_, placements_ ); }
private:
  unsigned int sequence_id_ = 0;
  std::vector<Placement> placements_;
  double entropy_ = -1.0;
};
