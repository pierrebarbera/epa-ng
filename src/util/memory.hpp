#pragma once

#include <string>

#include "core/raxml/Model.hpp"
#include "seq/MSA_Info.hpp"
#include "util/Options.hpp"
#include "util/logging.hpp"

class Options;

/*
  convert number of bytes into a formatted string of the highest sensible
  magnitude
 */
std::string format_byte_num( size_t size );

size_t memstring_to_byte( std::string s );

size_t get_max_memory();


/**
 * Class summarizing the estimated memory footprint of epa-ng, given the
 * specified inputs
 */
class Memory_Footprint {
  public:
  Memory_Footprint() = default;

  Memory_Footprint( MSA_Info const& ref_info,
                    MSA_Info const& qry_info,
                    raxml::Model const& model,
                    Options const& options );

  ~Memory_Footprint() = default;

  // getters
  size_t total() const
  {
    return partition_ + lookup_ + presample_ + refmsa_ + qsistream_ + allwork_;
  }

  size_t minimum() const
  {
    // size of partition assuming logn clv
    size_t const partition_logn = partition_ - clvbuffer_
        + ( logn_ * perclv_ );
    return partition_logn + presample_ + refmsa_ + qsistream_ + allwork_;
  }

  size_t clv() const { return perclv_; }
  size_t maximum_required_clvs() const { return maxnumclv_; }
  size_t logn_clvs() const { return logn_; }

  size_t lookup() const { return lookup_; }

  private:
  size_t partition_ = 0ul;
  size_t lookup_    = 0ul; // doubles as a signal that prescoring is enabled
  size_t presample_ = 0ul;
  size_t refmsa_    = 0ul;
  size_t qsistream_ = 0ul;
  size_t allwork_   = 0ul;

  size_t perclv_    = 0ul;
  size_t clvbuffer_ = 0ul;

  size_t maxnumclv_ = 0ul;
  size_t logn_      = 0ul;
};

/**
 * Class holding/setting the configuration necessary to govern how the memory
 * saver mode should allocate memory
 */
class Memory_Config {
  public:
  Memory_Config() = default;

  Memory_Config( std::string const& config_string ) { (void)config_string; }

  Memory_Config( Memory_Footprint const& footprint,
                 size_t const constraint )
  {
    auto const maxmem = get_max_memory();

    if( constraint > maxmem ) {
      LOG_WARN << "Specified memory limit of " << format_byte_num( constraint )
               << " exceeds the determined system wide maximum of "
               << format_byte_num( maxmem )
               << ". Continuing with the specified limit!";
    }

    auto const minmem = footprint.minimum();
    if( constraint < minmem ) {
      LOG_ERR << "Specified memory limit of " << format_byte_num( constraint )
              << " is below the minimum required value (for this input) of "
              << format_byte_num( minmem ) << ". Aborting!";
      std::exit( EXIT_FAILURE );
    }

    // how much above the minimum can we play with?
    // (this includes the logn + 2 clv slots)
    auto budget = constraint - minmem;

    // if we can afford it, use the preplacement lookup
    if( footprint.lookup() and ( footprint.lookup() < budget ) ) {
      budget -= footprint.lookup();
      preplace_lookup_enabled = true;
    }

    auto const per_clv = footprint.clv();

    // figure out how many more clv slots we can afford
    size_t extra_clv_slots = floor( static_cast< double >( budget ) / per_clv );

    // but have no more than the theoretical maximum
    clv_slots = std::min( footprint.logn_clvs() + extra_clv_slots,
                          footprint.maximum_required_clvs() );
  }

  ~Memory_Config() = default;

  // use clv_slots as indicator whether this config is valid
  operator bool() const { return clv_slots; }

  size_t concurrent_branches   = 4;
  bool preplace_lookup_enabled = true;
  size_t clv_slots             = 0u;
};

/**
 * High level option class administering the memsaver option
 */
class Memory_Saver {
  public:
  enum class Mode { kOff, kFull, kAuto, kCustom };

  Memory_Saver()  = default;

  Memory_Saver( Mode mode_arg,
                std::string maxmem_string,
                MSA_Info const& ref_info,
                MSA_Info const& qry_info,
                raxml::Model const& model,
                Options const& options )
      : mode( mode_arg )
      , footprint( ref_info, qry_info, model, options )
  {
    size_t const system_constraint = get_max_memory();
    size_t const user_constraint   = not maxmem_string.empty()
        ? memstring_to_byte( maxmem_string )
        : system_constraint;
    memory_constraint_ = std::min( system_constraint, user_constraint );

    init();
  }

  ~Memory_Saver() = default;

  Memory_Saver& operator=( Mode mode_arg )
  {
    // if the mode has changed, re-init
    if( mode != mode_arg ) {
      mode = mode_arg;
      init();
    } else {
      mode = mode_arg;
    }
    return *this;
  }

  Memory_Saver& operator=( bool enable )
  {
    if( enable ) {
      mode = Memory_Saver::Mode::kAuto;
    } else {
      mode = Memory_Saver::Mode::kOff;
    }
    return *this;
  }

  // no valid config = no memsave mode
  operator bool() const { return config; }

  Mode mode = Mode::kOff;
  Memory_Footprint footprint;
  Memory_Config config;

  private:
  void init()
  {
    switch( mode ) {
    case Mode::kCustom:
      std::runtime_error { "Custom memsave mode not implemented yet." };
      break;
    case Mode::kOff:
      break;
    case Mode::kAuto:
      // only create a valid memory config if we absolutely need it
      if( footprint.total() > memory_constraint_ * 0.95  ) {
        config = Memory_Config( footprint, memory_constraint_ );
      }
      break;
    case Mode::kFull:
      config = Memory_Config( footprint, footprint.minimum() );
      break;
    default:
      std::runtime_error { "Wrong mode!" };
    }
  }

  size_t memory_constraint_ = 0ul;
};
