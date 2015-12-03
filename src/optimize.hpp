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

void optimize(pll_utree_t * tree, pll_partition_t * partition, const Tree_Numbers &nums);

#endif
