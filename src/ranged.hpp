#pragma once

#include "Range.hpp"
#include "pllhead.hpp"

void update_partial_ranged(pll_partition_t * partition, pll_utree_t * tree, Range& range);
double compute_edge_logl_ranged(pll_partition_t * partition, pll_utree_t * tree, Range& range);
void fill_without(pll_partition_t * partition, unsigned int clv_index, Range& range, double value);
