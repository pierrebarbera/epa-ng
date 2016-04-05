#pragma once

#include <string>
#include <cstdio>
#include <stdexcept>

#include "pllhead.hpp"

void load_clv_from_binary(FILE* bin_file, long int offset, pll_partition_t * partition,
  const unsigned int clv_index)
{
  assert(bin_file);

  unsigned int attributes;
  auto err = pll_binary_clv_load(
    bin_file,
    0,
    partition,
    clv_index,
    &attributes,
    offset);

  if (err != PLL_SUCCESS)
    throw std::runtime_error{std::string("Loading CLV failed with pll_errno: ") + std::to_string(pll_errno)};
}

void load_tipchars_from_binary(FILE* bin_file, long int offset, pll_partition_t * partition,
  const unsigned int tipchar_index)
{
  (void)bin_file;
  (void)offset;
  (void)partition;
  (void)tipchar_index;
  // advance file pointer to correct tipchars

  // allocate space and read
  // auto ptr = (char*)calloc(tipchars_size, sizeof(char));
  // fread(ptr, sizeof(char), tipchars_size, bin_file);
  // return ptr;
}
