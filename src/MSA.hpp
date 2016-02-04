#pragma once

#include <vector>

#include "Sequence.hpp"

class MSA
{
public:
  typedef Sequence                                        value_type;
  typedef typename std::vector<Sequence>::iterator        iterator;
  typedef typename std::vector<Sequence>::const_iterator  const_iterator;

  MSA(const unsigned int num_sites);
  MSA();
  ~MSA();
  void move_sequences(iterator begin, iterator end);
  void append(const std::string& header, const std::string& sequence);
  void append(Sequence s);
  void erase(iterator begin, iterator end) {sequence_list_.erase(begin, end);};

  // getters
  unsigned int size() const {return sequence_list_.size();};
  unsigned int num_sites() const {return num_sites_;};
  const Sequence& operator[](const unsigned int i) const;

  // setters
  void num_sites(const unsigned int sites) {num_sites_ = sites;};

  //Iterator Compatability
  iterator begin() { return sequence_list_.begin(); };
  iterator end() { return sequence_list_.end(); };
  const_iterator begin() const { return sequence_list_.cbegin(); };
  const_iterator end() const { return sequence_list_.cend(); };
  const_iterator cbegin() { return sequence_list_.cbegin(); };
  const_iterator cend() { return sequence_list_.cend(); };

private:
  // Data Members
  unsigned int num_sites_;
  std::vector<Sequence> sequence_list_;
};
