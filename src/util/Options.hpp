#pragma once

#include <limits>
#include <string>

#include "util/memory.hpp"

/**
 * High level option class administering the memsaver option.
 *
 */
class Memsave_Option {
  public:
  enum class Mode { kOff, kFull, kAuto, kCustom };

  Memsave_Option()
      : Memsave_Option( Mode::kOff, "" )
  {
  }

  Memsave_Option( Mode mode_arg, std::string maxmem_string )
      : mode( mode_arg )
  {
    size_t const system_constraint = get_max_memory();
    size_t const user_constraint   = not maxmem_string.empty()
        ? memstring_to_byte( maxmem_string )
        : system_constraint;
    memory_constraint = std::min( system_constraint, user_constraint );
  }

  ~Memsave_Option() = default;

  Memsave_Option& operator=( Mode mode_arg )
  {
    mode = mode_arg;
    return *this;
  }

  Memsave_Option& operator=( bool enable )
  {
    if( enable ) {
      mode = Memsave_Option::Mode::kAuto;
    } else {
      mode = Memsave_Option::Mode::kOff;
    }
    return *this;
  }

  Mode mode                = Mode::kOff;
  size_t memory_constraint = 0ul;
};

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
  Numerical_Scaling scaling     = Numerical_Scaling::kAuto;
  bool preserve_rooting         = true;
  Memsave_Option memsave;
};
