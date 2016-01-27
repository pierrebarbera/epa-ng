#include "optimize.hpp"

#include <vector>
#include <cmath>
#include <cstdlib>
#include <cassert>
#include <stdexcept>
#include <limits>

#include "pll_util.hpp"
#include "constants.hpp"

using namespace std;

// local helpers
void traverse_update_partials(pll_utree_t * tree, pll_partition_t * partition,
  pll_utree_t ** travbuffer, double * branch_lengths, unsigned int * matrix_indices,
  pll_operation_t * operations);
double compute_negative_lnl_unrooted_ranged (void * p, double *x);
double optimize_triplet_lbfgsb_ranged (pll_partition_t * partition, pll_utree_t * tree,
  double factr, double pgtol, unsigned int begin, unsigned int span);
double compute_edge_loglikelihood_ranged(pll_partition_t * partition,
  unsigned int parent_clv_index, int parent_scaler_index, unsigned int child_clv_index,
  int child_scaler_index, unsigned int matrix_index, unsigned int freqs_index,
  unsigned int begin, unsigned int span);
void update_partial_ranged(pll_partition_t * partition, pll_operation_t * op,
  unsigned int begin, unsigned int span);

void update_partial_ranged(pll_partition_t * partition, pll_operation_t * op,
                            unsigned int begin, unsigned int span)
{
#ifdef __AVX
  auto attrib = PLL_ATTRIB_ARCH_AVX;
#else
  auto attrib = PLL_ATTRIB_ARCH_SSE;
#endif

  // make sure we actually skip the first <begin> entries ov the CLVs. Their size is
  // number of states * number of rate cats
  auto skip_clv = begin * partition->states * partition->rate_cats;
  // scalers are per site, so we just skip <begin>

  // check if scalers exist, and if so shift them accordingly
  unsigned int * parent_scaler = (op->parent_scaler_index == PLL_SCALE_BUFFER_NONE) ?
      nullptr : partition->scale_buffer[op->parent_scaler_index] + begin;
  unsigned int * child1_scaler = (op->child1_scaler_index == PLL_SCALE_BUFFER_NONE) ?
      nullptr : partition->scale_buffer[op->child1_scaler_index] + begin;
  unsigned int * child2_scaler = (op->child2_scaler_index == PLL_SCALE_BUFFER_NONE) ?
      nullptr : partition->scale_buffer[op->child2_scaler_index] + begin;

  pll_core_update_partial(partition->states,
                          span, // begin + span = first CLV entry not used in computation
                          partition->rate_cats,
                          partition->clv[op->parent_clv_index] + skip_clv,
                          parent_scaler,
                          partition->clv[op->child1_clv_index] + skip_clv,
                          partition->clv[op->child2_clv_index] + skip_clv,
                          partition->pmatrix[op->child1_matrix_index],
                          partition->pmatrix[op->child2_matrix_index],
                          child1_scaler,
                          child2_scaler,
                          attrib
                        );
}

/* Function to compute the edge likelihood directly on the CLVs passed
  instead of indirectly through the partition instance / indices
  for the purpose of ranged computation */
double compute_edge_loglikelihood_ranged(pll_partition_t * partition,
                                              unsigned int parent_clv_index,
                                              int parent_scaler_index,
                                              unsigned int child_clv_index,
                                              int child_scaler_index,
                                              unsigned int matrix_index,
                                              unsigned int freqs_index,
                                              unsigned int begin,
                                              unsigned int span)
{
  unsigned int n,i,j,k;
  double logl = 0;
  double terma, terma_r, termb;
  double site_lk, inv_site_lk;

  auto skip_clv = begin * partition->states * partition->rate_cats;
  auto pattern_weights = partition->pattern_weights + begin;
  auto invariant = partition->invariant + begin;

  const double * clvp = partition->clv[parent_clv_index] + skip_clv;
  const double * clvc = partition->clv[child_clv_index] + skip_clv;
  const double * freqs = NULL;
  const double * pmatrix = partition->pmatrix[matrix_index];
  double prop_invar = partition->prop_invar[freqs_index];
  unsigned int states = partition->states;
  unsigned int states_padded = partition->states_padded;
  double * weights = partition->rate_weights;
  unsigned int scale_factors;

  unsigned int * parent_scaler;
  unsigned int * child_scaler;

  if (child_scaler_index == PLL_SCALE_BUFFER_NONE)
    child_scaler = NULL;
  else
    child_scaler = partition->scale_buffer[child_scaler_index] + begin;

  if (parent_scaler_index == PLL_SCALE_BUFFER_NONE)
    parent_scaler = NULL;
  else
    parent_scaler = partition->scale_buffer[parent_scaler_index] + begin;

  for (n = 0; n < span; ++n)
  {
    pmatrix = partition->pmatrix[matrix_index];
    freqs = partition->frequencies[freqs_index];
    terma = 0;
    for (i = 0; i < partition->rate_cats; ++i)
    {
      terma_r = 0;
      for (j = 0; j < partition->states; ++j)
      {
        termb = 0;
        for (k = 0; k < partition->states; ++k)
          termb += pmatrix[k] * clvc[k];
        terma_r += clvp[j] * freqs[j] * termb;
        pmatrix += states;
      }
      terma += terma_r * weights[i];
      clvp += states;
      clvc += states;
      if (partition->mixture > 1)
        freqs += states_padded;
    }

    site_lk = terma;

    /* account for invariant sites */
    if (prop_invar > 0)
    {
      inv_site_lk = (invariant[n] == -1) ?
                        0 : freqs[invariant[n]];

      site_lk = site_lk * (1. - prop_invar) +
                inv_site_lk * prop_invar;
    }

    scale_factors = (parent_scaler) ? parent_scaler[n] : 0;
    scale_factors += (child_scaler) ? child_scaler[n] : 0;

    logl += log(site_lk) * pattern_weights[n];
    if (scale_factors)
      logl += scale_factors * log(PLL_SCALE_THRESHOLD);
  }

  return logl;
}

/* compute the negative lnL (score function for L-BFGS-B
* Original author:  Diego Darriba <Diego.Darriba@h-its.org>
* Used with permission.
* With modifications.
*/
double compute_negative_lnl_unrooted_ranged (void * p, double *x)
{
  lk_set * params = (lk_set *) p;
  pll_partition_t * partition = params->partition;
  double score;

  pll_update_prob_matrices (partition, 0,
                            params->matrix_indices,
                            x,
                            3);
  update_partial_ranged (partition, params->operation, params->begin, params->span);

  score = -1 * compute_edge_loglikelihood_ranged (
                partition,
                params->tree->clv_index,
                params->tree->scaler_index,
                params->tree->back->clv_index,
                params->tree->back->scaler_index,
                params->tree->pmatrix_index,
                0,
                params->begin,
                params->span);

  return score;
}

/* optimize the 3 adjacent branches of the node (tree) simultaneously
 * we assume that 'tree' is the virtual root and the CLVs of the nodes adjacent
 * to the optimized branches are valid.
 * This function returns the updated branch lengths in the tree structure and
 * the likelihood score.
 * Original author:  Diego Darriba <Diego.Darriba@h-its.org>
 * Used with permission.
 * With modifications.
 */
double optimize_triplet_lbfgsb_ranged (pll_partition_t * partition,
                                       pll_utree_t * tree,
                                       double factr,
                                       double pgtol,
                                       unsigned int begin,
                                       unsigned int span)
{

  /* L-BFGS-B parameters */
  unsigned int num_variables = 3;
  double score = 0;
  double x[3], lower_bounds[3], upper_bounds[3];
  int bound_type[3];

  /* set x vector */
  x[0] = tree->length;
  x[1] = tree->next->length;
  x[2] = tree->next->next->length;

  /* set boundaries */
  lower_bounds[0] = lower_bounds[1] = lower_bounds[2] = PLL_OPT_MIN_BRANCH_LEN;
  upper_bounds[0] = upper_bounds[1] = upper_bounds[2] = PLL_OPT_MAX_BRANCH_LEN;
  bound_type[0] = bound_type[1] = bound_type[2] = PLL_LBFGSB_BOUND_BOTH;

  /* set operation for updating the CLV on the virtual root edge */
  pll_operation_t op;
  op.child1_clv_index = tree->next->back->clv_index;
  op.child1_matrix_index = tree->next->back->pmatrix_index;
  op.child1_scaler_index = tree->next->back->scaler_index;
  op.child2_clv_index = tree->next->next->back->clv_index;
  op.child2_matrix_index = tree->next->next->back->pmatrix_index;
  op.child2_scaler_index = tree->next->next->back->scaler_index;
  op.parent_clv_index = tree->clv_index;
  op.parent_scaler_index = tree->scaler_index;

  /* set likelihood parameters */
  lk_set params;
  params.partition = partition;
  params.tree = tree;
  params.matrix_indices[0] = tree->pmatrix_index;
  params.matrix_indices[1] = tree->next->pmatrix_index;
  params.matrix_indices[2] = tree->next->next->pmatrix_index;
  params.operation = &op;
  params.begin = begin;
  params.span = span;

  /* optimize branches */
  score = pll_minimize_lbfgsb (x, lower_bounds, upper_bounds, bound_type,
                               num_variables, factr, pgtol, &params,
                               compute_negative_lnl_unrooted_ranged);

  /* set lengths back to the tree structure */
  tree->length = tree->back->length = x[0];
  tree->next->length = tree->next->back->length = x[1];
  tree->next->next->length = tree->next->next->back->length = x[2];

  return score;
}

double optimize_branch_triplet_ranged(pll_partition_t * partition, pll_utree_t * tree, Range range)
{
  // compute logl once to give us a logl starting point
  auto logl = -1* optimize_triplet_lbfgsb_ranged (partition, tree, 1e7, OPT_PARAM_EPSILON,
                                                range.begin, range.span);
  auto cur_logl = -numeric_limits<double>::infinity();

  while (fabs (cur_logl - logl) > OPT_EPSILON)
  {
    logl = cur_logl;
    cur_logl = -1* optimize_triplet_lbfgsb_ranged (partition, tree, 1e7, OPT_PARAM_EPSILON,
                                          range.begin, range.span);
  }

  return cur_logl;
}

void traverse_update_partials(pll_utree_t * tree, pll_partition_t * partition,
    pll_utree_t ** travbuffer, double * branch_lengths, unsigned int * matrix_indices,
    pll_operation_t * operations)
{
  unsigned int num_matrices, num_ops;
  /* perform a full traversal*/
  assert(tree->next != nullptr);
  unsigned int traversal_size;
  pll_utree_traverse(tree, cb_full_traversal, travbuffer, &traversal_size);

  /* given the computed traversal descriptor, generate the operations
     structure, and the corresponding probability matrix indices that
     may need recomputing */
  pll_utree_create_operations(travbuffer,
                              traversal_size,
                              branch_lengths,
                              matrix_indices,
                              operations,
                              &num_matrices,
                              &num_ops);

  pll_update_prob_matrices(partition,
                           0,             // use model 0
                           matrix_indices,// matrices to update
                           branch_lengths,
                           num_matrices); // how many should be updated

  /* use the operations array to compute all num_ops inner CLVs. Operations
     will be carried out sequentially starting from operation 0 towrds num_ops-1 */
  pll_update_partials(partition, operations, num_ops);

}

double optimize_branch_lengths(pll_utree_t * tree, pll_partition_t * partition, pll_optimize_options_t& params,
  pll_utree_t ** travbuffer, double cur_logl, double lnl_monitor, int* smoothings)
{
  if (!tree->next)
    tree = tree->back;

  traverse_update_partials(tree, partition, travbuffer, params.lk_params.branch_lengths,
    params.lk_params.matrix_indices, params.lk_params.operations);

  pll_errno = 0; // hotfix

  cur_logl = -1 * pll_optimize_branch_lengths_iterative(
      partition, tree, 0, 0,
      OPT_PARAM_EPSILON, *smoothings++, 1); // last param = 1 means branch lengths are iteratively
                                      // updated during the call
  if (cur_logl+1e-6 < lnl_monitor)
    throw runtime_error{string("cur_logl < lnl_monitor: ") + to_string(cur_logl) + string(" : ")
    + to_string(lnl_monitor)};

  // reupdate the indices as they may have changed
  params.lk_params.where.unrooted_t.parent_clv_index = tree->clv_index;
  params.lk_params.where.unrooted_t.parent_scaler_index =
      tree->scaler_index;
  params.lk_params.where.unrooted_t.child_clv_index = tree->back->clv_index;
  params.lk_params.where.unrooted_t.child_scaler_index =
      tree->back->scaler_index;
  params.lk_params.where.unrooted_t.edge_pmatrix_index =
      tree->pmatrix_index;

  return cur_logl;
}

void optimize(Model& model, pll_utree_t * tree, pll_partition_t * partition,
  const Tree_Numbers& nums, const bool opt_branches, const bool opt_model)
{
  if (!opt_branches && !opt_model)
    return;

  if (opt_branches)
    set_branch_length(tree, DEFAULT_BRANCH_LENGTH);

  auto symmetries = (&(model.symmetries())[0]);

  // sadly we explicitly need these buffers here and in the params structure
  vector<pll_utree_t*> travbuffer(nums.nodes);
  vector<double> branch_lengths(nums.branches);
  vector<unsigned int> matrix_indices(nums.branches);
  vector<pll_operation_t> operations(nums.nodes);

  traverse_update_partials(tree, partition, &travbuffer[0], &branch_lengths[0],
    &matrix_indices[0], &operations[0]);

  // compute logl once to give us a logl starting point
  auto logl = pll_compute_edge_loglikelihood (partition, tree->clv_index,
                                                tree->scaler_index,
                                                tree->back->clv_index,
                                                tree->back->scaler_index,
                                                tree->pmatrix_index, 0);

  double cur_logl = -numeric_limits<double>::infinity();
  int smoothings = 1;
  double lnl_monitor = logl;

  // set up high level options structure
  pll_optimize_options_t params;
  params.lk_params.partition = partition;
  params.lk_params.operations = &operations[0];
  params.lk_params.branch_lengths = &branch_lengths[0];
  params.lk_params.matrix_indices = &matrix_indices[0];
  params.lk_params.alpha_value = model.alpha();
  params.lk_params.freqs_index = 0;
  params.lk_params.rooted = 0;
  params.lk_params.where.unrooted_t.parent_clv_index = tree->clv_index;
  params.lk_params.where.unrooted_t.parent_scaler_index = tree->scaler_index;
  params.lk_params.where.unrooted_t.child_clv_index = tree->back->clv_index;
  params.lk_params.where.unrooted_t.child_scaler_index =
    tree->back->scaler_index;
  params.lk_params.where.unrooted_t.edge_pmatrix_index = tree->pmatrix_index;

  /* optimization parameters */
  params.params_index = 0;
  params.mixture_index = 0;
  params.subst_params_symmetries = symmetries;
  params.factr = 1e7;
  params.pgtol = OPT_PARAM_EPSILON;

  while (fabs (cur_logl - logl) > OPT_EPSILON)
  {
    logl = cur_logl;

    if (opt_model)
    {
      params.which_parameters = PLL_PARAMETER_ALPHA;
      pll_optimize_parameters_brent(&params);

      params.which_parameters = PLL_PARAMETER_SUBST_RATES;
      pll_optimize_parameters_lbfgsb(&params);

      // params.which_parameters = PLL_PARAMETER_FREQUENCIES;
      // pll_optimize_parameters_lbfgsb(&params);

      params.which_parameters = PLL_PARAMETER_PINV;
      cur_logl = -1 * pll_optimize_parameters_brent(&params);
    }

    if (opt_branches)
      cur_logl = optimize_branch_lengths(tree,
          partition, params, &travbuffer[0], cur_logl, lnl_monitor, &smoothings);
  }

  if (opt_model)
  {
    // update epa model object as well
    model.alpha(params.lk_params.alpha_value);
    model.substitution_rates(partition->subst_params[0], 6);
    model.base_frequencies(partition->frequencies[params.params_index], partition->states);
  }
}

void compute_and_set_empirical_frequencies(pll_partition_t * partition, Model& model)
{
  double * empirical_freqs = pll_compute_empirical_frequencies (partition);
  pll_set_frequencies (partition, 0, 0, empirical_freqs);
  model.base_frequencies(partition->frequencies[0], partition->states);
  free (empirical_freqs);
}
