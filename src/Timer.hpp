#pragma once

#include <chrono>
#include <vector>
#include <algorithm>

class Timer {
  typedef typename std::chrono::milliseconds duration;
  typedef typename std::chrono::high_resolution_clock clock;
public:
  Timer() = default;
  ~Timer() = default;

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
private:
  clock::time_point start_;
  std::vector<duration> ts_;
};
