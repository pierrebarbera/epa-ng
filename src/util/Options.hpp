#pragma once

#include <limits>
#include <string>

class MemoryConfig {
  public:
  MemoryConfig()  = default;

  MemoryConfig( std::string const& config_string )
  {
    // build from a config string
    (void)config_string;
  }

  ~MemoryConfig() = default;

  size_t concurrent_branches   = 10;
  bool preplace_lookup_enabled = true;
};

class Options {

  public:
  enum class NumericalScaling { kOn, kOff, kAuto };

  class MemorySaver {
  public:
    enum class Mode { kOff, kFull, kAuto, kCustom };

    MemorySaver() = default;
    ~MemorySaver() = default;

    MemorySaver& operator=( Mode mode_ )
    {
      mode = mode_;
      return *this;
    }
    MemorySaver& operator=( bool enable )
    {
      if( enable ) {
        mode = Options::MemorySaver::Mode::kFull;
      } else {
        mode = Options::MemorySaver::Mode::kOff;
      }
      return *this;
    }

    operator bool() const { return mode != Mode::kOff; }
    // operator int() const { return static_cast<int>(bool()); }

    Mode mode = Mode::kOff;
  };

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
  unsigned int chunk_size       = 5000;
  unsigned int num_threads      = 0;
  bool repeats                  = false;
  bool premasking               = true;
  bool baseball                 = false;
  std::string tmp_dir;
  unsigned int precision        = 10;
  NumericalScaling scaling      = NumericalScaling::kAuto;
  bool preserve_rooting         = true;
  MemorySaver memsave;

  MemoryConfig memory_config;
};
