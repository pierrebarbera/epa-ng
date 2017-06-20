#pragma once

#include <limits>

class Options {

public:
  Options()  = default;
  ~Options() = default;

  bool prescoring               = false;
  bool opt_model                = false;
  bool opt_branches             = false;
  bool sliding_blo              = true;
  double support_threshold      = 0.9999;
  bool acc_threshold            = true;
  unsigned int filter_min       = 1;
  unsigned int filter_max       = std::numeric_limits<unsigned int>::max();
  bool prescoring_by_percentage = false;
  double prescoring_threshold   = 0.95;
  bool ranged                   = false;
  bool dump_binary_mode         = false;
  bool load_binary_mode         = false;
  unsigned int chunk_size       = 1000;
  unsigned int num_threads      = 0;
  bool repeats                  = true;
};
