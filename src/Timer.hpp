#pragma once

#include <chrono>
#include <vector>
#include <algorithm>
#include <cereal/types/vector.hpp>
#include <cereal/types/chrono.hpp>

class Timer {
  typedef typename std::chrono::milliseconds              duration;
  typedef typename std::chrono::high_resolution_clock     clock;
  typedef typename std::vector<duration>::iterator        iterator;
  typedef typename std::vector<duration>::const_iterator  const_iterator;
public:
  Timer() = default;
  ~Timer() = default;

  // Iterator Compatibility
  iterator begin() { return ts_.begin(); }
  iterator end() { return ts_.end(); }
  const_iterator begin() const { return ts_.cbegin(); }
  const_iterator end() const { return ts_.cend(); }
  const_iterator cbegin() { return ts_.cbegin(); }
  const_iterator cend() { return ts_.cend(); }

  void insert(iterator position, const_iterator first, const_iterator last)
  { ts_.insert(position, first, last); };

  void start() {start_ = clock::now();}
  void stop()
  {
    auto end = clock::now();
    auto runtime = std::chrono::duration_cast<duration>(end - start_);
    ts_.push_back(runtime);
  }

  double average()
  {
    duration sum(0);
    duration size(ts_.size());
    for (auto p : ts_)
      sum += p;
    return sum/size;
  }

  void clear() {ts_.clear();}

  // serialization
  template <class Archive>
  void save(Archive & ar) const { ar( ts_ ); }

  template <class Archive>
  void load(Archive & ar) { ar( ts_ ); }
private:
  clock::time_point start_;
  std::vector<duration> ts_;
};
