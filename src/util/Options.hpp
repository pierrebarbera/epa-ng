#pragma once

#include <limits>
#include <string>

#include "util/memory.hpp"

class Options {

  public:
  enum class Numerical_Scaling { kOn, kOff, kAuto };

  Options()  = default;
  ~Options() = default;

  bool prescoring               = true;
  bool opt_model                = false;
  bool opt_branches             = false;
  bool sliding_blo              = true;
  double support_threshold      = 0.01;
  bool acc_threshold            = false;
  unsigned int filter_min       = 1;
  unsigned int filter_max       = 7;
  bool prescoring_by_percentage = false;
  double prescoring_threshold   = 0.99999;
  bool ranged                   = false;
  bool dump_binary_mode         = false;
  bool load_binary_mode         = false;
  size_t chunk_size             = 5000;
  unsigned int num_threads      = 0;
  bool repeats                  = false;
  bool premasking               = true;
  bool baseball                 = false;
  std::string tmp_dir;
  unsigned int precision        = 10;
  Numerical_Scaling scaling      = Numerical_Scaling::kAuto;
  bool preserve_rooting         = true;
  Memory_Saver memsave;
};
