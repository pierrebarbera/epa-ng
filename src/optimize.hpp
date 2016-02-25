#pragma once

#include "pllhead.hpp"
#include "Model.hpp"
#include "Tree_Numbers.hpp"
#include "Range.hpp"

constexpr double OPT_EPSILON = 1.0;
constexpr double OPT_PARAM_EPSILON = 1e-4;
constexpr double OPT_BRANCH_EPSILON = 1e-4;
// #define OPT_BRANCH_EPSILON  -log(1e-4)
constexpr double OPT_FACTR = 1e7;

typedef struct
{
  pll_utree_t * tree;
  pll_partition_t * partition;
  unsigned int matrix_indices[3];
  pll_operation_t * operation;
  unsigned int begin;
  unsigned int span;
}lk_set;

// interface
void optimize(Model& model, pll_utree_t * tree, pll_partition_t * partition,
  const Tree_Numbers& nums, const bool opt_branches, const bool opt_model);
double optimize_branch_triplet_ranged(pll_partition_t * partition, pll_utree_t * tree, Range range);
void compute_and_set_empirical_frequencies(pll_partition_t * partition, Model& model);
double optimize_branch_triplet_newton(pll_partition_t * partition, pll_utree_t * tree);
