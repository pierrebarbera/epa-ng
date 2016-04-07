#pragma once

#include "pllhead.hpp"

void tiny_partition_destroy(pll_partition_t * partition);
pll_utree_t * make_tiny_tree_structure(const pll_utree_t * old_proximal, const pll_utree_t * old_distal,
  const bool tip_tip_case);
pll_partition_t * make_tiny_partition(const pll_partition_t * old_partition, const pll_utree_t * tree,
  const pll_utree_t * old_proximal, const pll_utree_t * old_distal, const bool tip_tip_case);
