#pragma once

#include <string>
#include <cstdio>

#include "pllhead.hpp"

class Binary {
public:
  Binary(const std::string& bin_file_path);
  Binary() : bin_fptr_(nullptr) { }
  ~Binary() { if(bin_fptr_) fclose(bin_fptr_); }

  // access functions
  void load_clv(pll_partition_t * partition, const unsigned int clv_index);
  void load_tipchars(pll_partition_t * partition, const unsigned int tipchars_index);
  void load_scaler(pll_partition_t * partition, const unsigned int scaler_index);
  pll_partition_t* load_partition();
  pll_utree_t* load_utree();
private:
  FILE* bin_fptr_;
};

class Tree;

void dump_to_binary(Tree& tree, const std::string& file);
