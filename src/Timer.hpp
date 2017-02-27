#pragma once

#include <chrono>
#include <vector>
#include <algorithm>
#include <cereal/types/vector.hpp>
#include <cereal/types/chrono.hpp>

class Timer {
public:
  // Typedefs
  using duration        = std::chrono::milliseconds;
  using clock           = std::chrono::high_resolution_clock;
  using iterator        = std::vector<duration>::iterator;
  using const_iterator  = std::vector<duration>::const_iterator;
  
  // Constructors/Destructors
  Timer(std::vector<double> init_list) 
  {
    for (auto elem : init_list)
    {
      duration fp_ms(static_cast<unsigned int>(elem*1000));
      ts_.push_back(fp_ms);
    }
  }
  Timer()   = default;
  ~Timer()  = default;

  // Iterator Compatibility
  iterator begin() { return ts_.begin(); }
  iterator end() { return ts_.end(); }
  const_iterator begin() const { return ts_.cbegin(); }
  const_iterator end() const { return ts_.cend(); }
  const_iterator cbegin() { return ts_.cbegin(); }
  const_iterator cend() { return ts_.cend(); }

  // Methods
  void insert(iterator position, const_iterator first, const_iterator last)
  { ts_.insert(position, first, last); };

  void start() 
  {
    start_ = clock::now();
  }

  void pause()
  {
    pause_start_ = clock::now();
  }

  void resume()
  {
    auto end = clock::now();
    auto pause_time = std::chrono::duration_cast<duration>(end - pause_start_);
    pauses_.push_back(pause_time);
  }

  void stop()
  {
    auto end = clock::now();
    duration pause_total(0);
    for (auto p : pauses_)
      pause_total += p;

    auto runtime = std::chrono::duration_cast<duration>(end - start_) - pause_total;
    ts_.push_back(runtime);
    pauses_.clear();
  }

  double sum()
  {
    duration sum(0);
    for (auto p : ts_)
      sum += p;
    return sum.count();
  }

  double average()
  {
    return this->sum()/ts_.size();
  }

  void clear() {ts_.clear();}

  // Serialization
  template <class Archive>
  void save(Archive & ar) const { ar( ts_ ); }

  template <class Archive>
  void load(Archive & ar) { ar( ts_ ); }
private:
  clock::time_point start_;
  std::vector<duration> ts_;
  clock::time_point pause_start_;
  std::vector<duration> pauses_;
};
