#pragma once

#include <algorithm>
#include <unordered_map>
#include <cereal/types/unordered_map.hpp>

#include "Sample.hpp"

/**
 * work[branch_id] = {seq_id_1. seq_id_2, ...}
 */
class Work
{
public:
  typedef unsigned int                              key_type;
  typedef std::vector<unsigned int>                 value_type;
  typedef std::unordered_map<key_type, value_type>  container_type;
  typedef typename container_type::iterator         iterator;
  typedef typename container_type::const_iterator   const_iterator;

  /**
   * Create work object from a Sample: all entries are seen as placements to be recomputed
   */
  Work(Sample& sample)
  {
    for (auto & pq : sample)
    {
      auto seq_id = pq.sequence_id();
      for (auto & placement : pq)
      {
        auto branch_id = placement.branch_id();
        work_set_[branch_id].push_back(seq_id);
      }
    }
  }

  /**
   * Create a work object covering all sequences in [seq_range.first, seq_range.second)
   * for every branch ID in [branch_range.first, branch_range.second).
   */
  Work(std::pair<size_t, size_t>&& branch_range, std::pair<size_t, size_t>&& seq_range)
  {
    for (size_t branch_id = branch_range.first; branch_id < branch_range.second; ++branch_id)
      for (size_t seq_id = seq_range.first; seq_id < seq_range.second; ++seq_id)
        work_set_[branch_id].push_back(seq_id);      
  }

  Work() = default;       

  ~Work() = default;

  void clear() {work_set_.clear();}
  size_t size() const 
  {
    return std::accumulate(work_set_.begin(), work_set_.end(), 0,
      [](size_t a, std::pair<key_type, value_type> b){return a + b.second.size();}
      );
  }

  // template< class InputIt >
  // void insert(InputIt first, InputIt last) {work_set_.insert(first, last);}

  // Iterator Compatibility
  iterator begin() { return work_set_.begin(); }
  iterator end() { return work_set_.end(); }
  const_iterator begin() const { return work_set_.cbegin(); }
  const_iterator end() const { return work_set_.cend(); }
  const_iterator cbegin() { return work_set_.cbegin(); }
  const_iterator cend() { return work_set_.cend(); }

  // Operator overloads
  value_type& operator[] (const unsigned int index) { return work_set_[index]; }
  
  // serialization
  template <class Archive>
  void save(Archive & ar) const { ar( work_set_ ); }

  template <class Archive>
  void load(Archive & ar) { ar( work_set_ ); }

private:
  container_type work_set_;
};
