#ifndef EPA_OPTIMIZE_H_
#define EPA_OPTIMIZE_H_

#include "pllhead.hpp"
#include "Model.hpp"
#include "Tree_Numbers.hpp"
#include "Range.hpp"

constexpr double OPT_EPSILON = 1.0;
constexpr double OPT_PARAM_EPSILON = 1e-2;

/* if set, the parameters are no longer optimized when
 * they do not improve the likelihood at one iteration */
#define CHECK_LOCAL_CONVERGENCE 1

typedef struct
{
  pll_utree_t * tree;
  pll_partition_t * partition;
  unsigned int matrix_indices[3];
  pll_operation_t * operation;
}lk_set;

// local helpers
void traverse_update_partials(pll_utree_t * tree, pll_partition_t * partition,
    pll_utree_t ** travbuffer, double * branch_lengths, unsigned int * matrix_indices,
    pll_operation_t * operations);

// interface
double optimize_branch_lengths(pll_utree_t * tree, pll_partition_t * partition, pll_optimize_options_t& params,
  pll_utree_t ** travbuffer, double cur_log, double lnl_monitor, int* smoothings);
void optimize(Model& model, pll_utree_t * tree, pll_partition_t * partition,
  const Tree_Numbers& nums, const bool opt_branches, const bool opt_model);
double optimize_branch_triplet(pll_partition_t * partition, pll_utree_t * tree);
void compute_and_set_empirical_frequencies(pll_partition_t * partition, Model& model);

#endif
