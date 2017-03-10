#pragma once

#include <mutex>
#include <vector>
#include <map>
#include <utility>
#include <limits>

#include "Matrix.hpp"

static constexpr size_t NT_MAP_SIZE = 16;
static constexpr size_t AA_MAP_SIZE = 25;
static constexpr unsigned char NT_MAP[16] = 
  {'A', 'C', 'G', 'T', '-', 'Y', 'R', 'W', 'S', 'K', 'M', 'D', 'V', 'H', 'B', 'X'};
static constexpr unsigned char AA_MAP[25] = 
  { 'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 
    'T', 'V', 'W', 'Y', '-', 'X', 'B', 'Z', 'J'};
static constexpr size_t INVALID = std::numeric_limits<size_t>::max();

class Lookup_Store
{
public:
  using lookup_type = Matrix<double>;

  Lookup_Store(const size_t num_branches, const size_t num_states) 
    : branch_(num_branches)
    , store_(num_branches)
    , char_map_size_((num_states == 4) ? NT_MAP_SIZE : AA_MAP_SIZE)
    , char_map_((num_states == 4) ? NT_MAP : AA_MAP)
  {
    for (size_t i = 0; i < 128; ++i) {
      char_to_posish_[i] = INVALID;
    }
    // build reverse map from char to ID in the charmap
    for (size_t i = 0; i < char_map_size_; ++i) {
      char_to_posish_[char_map_[i]] = i;
    }
  }

  Lookup_Store() = delete;
  ~Lookup_Store() = default;

  void init_branch(const size_t branch_id, std::vector<std::vector<double>> precomps)
  {
    store_[branch_id] = Matrix<double>(precomps[0].size(), char_map_size_);

    for(size_t ch = 0; ch < precomps.size(); ++ch) {
      for(size_t site = 0; site < precomps[ch].size(); ++site) {
        store_[branch_id](site, ch) = precomps[ch][site];
      }
    }
  }

  std::mutex& get_mutex(const size_t branch_id)
  {
    return branch_[branch_id];
  }

  bool has_branch(const size_t branch_id) 
  {
    return store_[branch_id].size() != 0; 
  }

  lookup_type& operator[](const size_t branch_id)
  {
    return store_[branch_id];
  }

  unsigned char char_map(const size_t i)
  {
    if (i >= char_map_size_) {
      throw std::runtime_error{
        std::string("char_map access out of bounds! i =") + std::to_string(i)
      };
    }

    return char_map_[i];
  }

  size_t char_map_size()
  {
    return char_map_size_;
  }

  size_t char_position(unsigned char c) const
  {
    auto pos = char_to_posish_[c];

    if (pos == INVALID) {
      throw std::runtime_error{
        std::string("char is invalid! char =") + std::to_string(c)
      };
    }

    return pos;
  }

  double sum_precomputed_sitelk(const size_t branch_id, const std::string& seq) const
  {
    #ifndef NDEBUG
    assert(seq.length() == store_[branch_id].rows());
    #endif
    
    double sum = 0;
    for (size_t i = 0; i < seq.length(); ++i)
    {
      auto position = this->char_position(seq[i]);
      sum += store_[branch_id](i, position);
    }
    return sum;
  }

private:
  std::vector<std::mutex> branch_;
  std::vector<lookup_type> store_;
  const size_t char_map_size_;
  const unsigned char * char_map_;
  std::array<size_t, 128> char_to_posish_;
};
