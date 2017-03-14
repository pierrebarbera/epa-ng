#pragma once

#include <vector>
#include <utility>

#include "Sequence.hpp"

class MSA
{
public:
  using container_type  = std::vector<Sequence>;
  using iterator        = container_type::iterator;
  using const_iterator  = container_type::const_iterator;

  MSA(const size_t num_sites) : num_sites_(num_sites) {};
  MSA() : num_sites_(0) {};
  ~MSA() = default;

  void move_sequences(iterator begin, iterator end);
  void append(const std::string& header, const std::string& sequence);
  // void append(Sequence s);
  void erase(iterator begin, iterator end) {sequence_list_.erase(begin, end);}
  void clear() {sequence_list_.clear();}

  // getters
  size_t size() const {return sequence_list_.size();}
  size_t num_sites() const {return num_sites_;}
  const Sequence& operator[](const size_t i) const;

  // setters
  void num_sites(const size_t sites) {num_sites_ = sites;}

  //Iterator Compatability
  iterator begin() { return sequence_list_.begin(); };
  iterator end() { return sequence_list_.end(); };
  const_iterator begin() const { return sequence_list_.cbegin(); };
  const_iterator end() const { return sequence_list_.cend(); };
  const_iterator cbegin() const { return sequence_list_.cbegin(); };
  const_iterator cend() const { return sequence_list_.cend(); };

private:
  size_t num_sites_;
  container_type sequence_list_;
};
