#pragma once

#include <string>

#include "core/raxml/Model.hpp"
#include "core/pll/pll_util.hpp"
#include "core/pll/error.hpp"
#include "seq/MSA_Info.hpp"
#include "util/logging.hpp"


/*
  convert number of bytes into a formatted string of the highest sensible
  magnitude
 */
std::string format_byte_num( size_t size );

size_t memstring_to_byte( std::string s );

size_t get_max_memory();


#include "util/Options.hpp"
class Options;
class Memsave_Option;

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

  size_t minimum( size_t const concurrent_branches = 0 ) const
  {
    assert( partition_ );
    // size of partition assuming logn clv
    assert( partition_ > clvbuffer_ );
    size_t const partition_logn = partition_ - clvbuffer_ + ( logn_ * perclv_ );

    return partition_logn + presample_ + refmsa_ + qsistream_ + allwork_
        + memsaver_ + (concurrent_branches * perttdeep_);
  }

  size_t clv() const { return perclv_; }
  size_t maximum_required_clvs() const { return maxnumclv_; }
  size_t logn_clvs() const { return logn_; }
  size_t lookup() const { return lookup_; }

  size_t per_deepcopy_tiny_trees() const { return perttdeep_; }

  operator bool() const { return partition_; }

  private:
  size_t partition_ = 0ul;
  size_t lookup_    = 0ul; // doubles as a signal that prescoring is enabled
  size_t presample_ = 0ul;
  size_t refmsa_    = 0ul;
  size_t qsistream_ = 0ul;
  size_t allwork_   = 0ul;
  size_t tinytrees_ = 0ul;
  size_t memsaver_  = 0ul;

  size_t perclv_    = 0ul;
  size_t clvbuffer_ = 0ul;
  size_t perttdeep_ = 0ul;

  size_t maxnumclv_ = 0ul;
  size_t logn_      = 0ul;
};

/**
 * Extra structures needed when using the pll partition memory saving mode
 */
class Logn_Structures {
  public:
  Logn_Structures() = default;
  Logn_Structures( pll_utree_t* tree )
      : subtree_sizes_( pll_utree_get_subtree_sizes( tree ), free )
      , traversal_{ tree->edge_count, nullptr }
  {
    handle_pll_failure( not subtree_sizes_,
                        "pll_utree_get_subtree_sizes failed." );

    // get the traversal, hopefully one that minimizes overall recomputations
    utree_query_branches( tree, &traversal_[ 0 ] );
  }

  ~Logn_Structures() = default;

  Logn_Structures( Logn_Structures const& other ) = delete;
  Logn_Structures( Logn_Structures&& other )      = default;

  Logn_Structures& operator=( Logn_Structures const& other ) = delete;
  Logn_Structures& operator=( Logn_Structures&& other ) = default;

  operator bool() const { return subtree_sizes_ and not traversal_.empty(); }

  unsigned int const* subtree_sizes() const { return subtree_sizes_.get(); }
  std::vector< pll_unode_t* > const& traversal() const { return traversal_; }
  pll_unode_t const* traversal( size_t i ) const { return traversal_[ i ]; }
  pll_unode_t* traversal( size_t i ) { return traversal_[ i ]; }
  void reverse_traversal()
  {
    std::reverse( std::begin( traversal_ ), std::end( traversal_ ) );
  };

  private:
  std::unique_ptr< unsigned int, decltype( free )* > subtree_sizes_{ nullptr,
                                                                     free };
  std::vector< pll_unode_t* > traversal_;
};

/**
 * Class holding/setting the configuration necessary to govern how the memory
 * saver mode should allocate memory
 */
class Memory_Config {
  public:
  Memory_Config() = default;

  Memory_Config(  Memsave_Option const& memsave_opt,
                  Memory_Footprint const& footprint,
                  pll_utree_t* tree );

  ~Memory_Config() = default;

  Memory_Config( Memory_Config const& other ) = delete;
  Memory_Config( Memory_Config&& other )      = default;

  Memory_Config& operator=( Memory_Config const& other ) = delete;
  Memory_Config& operator=( Memory_Config&& other ) = default;

  // use clv_slots and validity of the structs as indicator whether this config
  // is valid
  operator bool() const { return clv_slots and structs; }

  size_t concurrent_branches   = 4;
  bool preplace_lookup_enabled = true;
  size_t clv_slots             = 0u;

  Logn_Structures structs;
  private:
  void init( size_t const constraint,
             Memory_Footprint const& footprint,
             pll_utree_t* tree );
};
