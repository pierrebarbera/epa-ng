#pragma once

#include <string>
#include <cstdio>
#include <vector>
#include <memory>

#include "pllhead.hpp"

typedef std::unique_ptr<FILE, int(*)(FILE*)> unique_fptr;
// custom deleter
int safe_fclose(FILE* fptr);

class Binary {
public:
  Binary(const std::string& bin_file_path);
  Binary() : bin_fptr_(nullptr, safe_fclose) { }
  Binary(Binary && other) = default;
  ~Binary() = default;

  Binary& operator=(Binary && other) = default;

  // access functions
  void load_clv(pll_partition_t * partition, const unsigned int clv_index);
  void load_tipchars(pll_partition_t * partition, const unsigned int tipchars_index);
  void load_scaler(pll_partition_t * partition, const unsigned int scaler_index);
  pll_partition_t* load_partition();
  pll_utree_t* load_utree();
private:
  unique_fptr bin_fptr_;
  std::vector<pll_block_map_t> map_;
};

class Tree;

void dump_to_binary(Tree& tree, const std::string& file);
