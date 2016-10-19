#pragma once

#include <string>
#include <vector>
#include <algorithm>
#include <cereal/types/vector.hpp>

#include "PQuery.hpp"

class Sample {
public:
  typedef PQuery                                            value_type;
  typedef typename std::vector<value_type>::iterator        iterator;
  typedef typename std::vector<value_type>::const_iterator  const_iterator;

  Sample() = default;
  Sample(unsigned int size) : pquerys_(size) {};
  Sample(const std::string newick) : newick_(newick) {};
  ~Sample() = default;

  // member access
  PQuery& back() { return pquerys_.back(); }
  unsigned int size() const { return pquerys_.size(); }
  const std::string& newick() const { return newick_; }
  void clear() { pquerys_.clear(); }
  void push_back(PQuery&& pq) { pquerys_.push_back(pq); }
  void push_back(PQuery& pq) { pquerys_.push_back(pq); }
  void erase(iterator begin, iterator end) { pquerys_.erase(begin, end); }

  // needs to be in the header
  template <typename ...Args>
  void emplace_back(Args && ...args) { pquerys_.emplace_back(std::forward<Args>(args)...); }
  template <class InputIt>
  void insert(InputIt first, InputIt last) {pquerys_.insert(pquerys_.end(), first, last);}

  template <typename ...Args>
  void add_placement(unsigned int seq_id, Args&& ...args)
  {
    // if seq_id in pquerys_
    auto iter = std::end(pquerys_);
    if ((iter = std::find(std::begin(pquerys_), std::end(pquerys_), PQuery(seq_id))) 
      != std::end(pquerys_))
    {
      iter->emplace_back(std::forward<Args>(args)...);
    }
    else
    {
      pquerys_.push_back(seq_id);
      pquerys_.back().emplace_back(std::forward<Args>(args)...);
    }
  }

  // Iterator Compatibility
  iterator begin() { return pquerys_.begin(); }
  iterator end() { return pquerys_.end(); }
  const_iterator begin() const { return pquerys_.cbegin(); }
  const_iterator end() const { return pquerys_.cend(); }
  const_iterator cbegin() { return pquerys_.cbegin(); }
  const_iterator cend() { return pquerys_.cend(); }

  // Operator overloads
  PQuery& operator[] (const unsigned int index) { return pquerys_[index]; }

  // serialization
  template <class Archive>
  void save(Archive & ar) const { ar( pquerys_ ); }

  template <class Archive>
  void load(Archive & ar) { ar( pquerys_ ); }

private:
  std::vector<PQuery> pquerys_;
  std::string newick_;
};
