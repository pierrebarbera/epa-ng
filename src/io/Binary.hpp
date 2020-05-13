#pragma once

#include <cstdio>
#include <memory>
#include <mutex>
#include <string>
#include <vector>

#include "core/pll/pllhead.hpp"

// custom deleter
int safe_fclose( FILE* fptr );

class Binary {
  public:
  using file_ptr_type = std::unique_ptr< FILE, int ( * )( FILE* ) >;

  explicit Binary( std::string const& bin_file_path );
  Binary()
      : bin_fptr_( nullptr, safe_fclose )
  {
  }
  Binary( Binary&& other );
  ~Binary() = default;

  Binary& operator=( Binary&& other );

  // access functions
  void load_clv( pll_partition_t* partition, unsigned int const clv_index );
  void load_tipchars( pll_partition_t* partition, unsigned int const tipchars_index );
  void load_scaler( pll_partition_t* partition, unsigned int const scaler_index );
  pll_partition_t* load_partition();
  pll_utree_t* load_utree( unsigned int const num_tips );

  private:
  std::mutex file_mutex_;
  file_ptr_type bin_fptr_;
  std::vector< pll_block_map_t > map_;
};

class Tree;

void dump_to_binary( Tree& tree, std::string const& file );
