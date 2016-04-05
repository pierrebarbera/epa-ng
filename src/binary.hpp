#pragma once

#include <string>
#include <cstdio>

#include "pllhead.hpp"

class Binary {
public:
  Binary(const std::string& bin_file_path);
  Binary() : bin_fptr_(nullptr) { }
  ~Binary() { fclose(bin_fptr_); }

  // access functions
  void load_clv(pll_partition_t * partition, const unsigned int clv_index);
  void load_tipchars(pll_partition_t * partition, const unsigned int tipchars_index);
  pll_partition_t* load_partition();
private:
  FILE* bin_fptr_;
  long int tipchars_offset_;
  long int clv_offset_;
  long int partition_offset_;
};
