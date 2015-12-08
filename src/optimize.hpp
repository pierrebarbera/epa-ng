#ifndef EPA_OPTIMIZE_H_
#define EPA_OPTIMIZE_H_

#include "pllhead.hpp"
#include "Model.hpp"
#include "Tree_Numbers.hpp"

#define OPT_EPSILON       1e-1
#define OPT_PARAM_EPSILON 1e-1

/* if set, the parameters are no longer optimized when
 * they do not improve the likelihood at one iteration */
#define CHECK_LOCAL_CONVERGENCE 1

// local helpers
void traverse_update_partials(pll_utree_t * tree, pll_partition_t * partition,
    pll_utree_t ** travbuffer, double * branch_lengths, unsigned int * matrix_indices,
    pll_operation_t * operations);

// interface
void optimize_branch_lengths(pll_utree_t * tree, pll_partition_t * partition, const Tree_Numbers &nums);
void optimize_model_params(Model& model, pll_utree_t * tree, pll_partition_t * partition,
  const Tree_Numbers& nums);

#endif
