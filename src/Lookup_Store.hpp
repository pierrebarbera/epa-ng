#pragma once

#include <mutex>
#include <vector>
#include <utility>

class Lookup_Store
{
public:
  typedef typename std::vector<std::vector<double>>     lookup_t;

  Lookup_Store(size_t num_branches) : branch_(num_branches), store_(num_branches) { }
  Lookup_Store() = default;
  ~Lookup_Store() = default;

  std::mutex& get_mutex(const size_t branch_id)
  {
    return branch_[branch_id];
  }

  bool has_branch(const size_t branch_id) 
  {
    return store_[branch_id].size() != 0; 
  }

  lookup_t& operator[](const size_t branch_id)
  {
    return store_[branch_id];
  }

private:
  std::vector<std::mutex> branch_;
  std::vector<lookup_t> store_;

};