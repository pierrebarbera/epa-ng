#include "optimize.hpp"

#include <vector>
#include <cmath>
#include <cstdlib>
#include <cassert>
#include <stdexcept>
#include <limits>

#include "pll_util.hpp"
#include "constants.hpp"
#include "logging.hpp"
#include "ranged.hpp"

using namespace std;

// TODO change to compute triplet edge logl? also make ranged
static void update_matrices_partial_triplet(pll_partition_t * partition, pll_utree_t * tree, Range& range)
{
  auto parent = tree;
  auto child1 = parent->next->back;
  auto child2 = parent->next->next->back;
  // update p-matrices
  double branch_lengths[3] = {parent->length, child1->length, child2->length};
  unsigned int matrix_indices[3] = {parent->pmatrix_index, child1->pmatrix_index, child2->pmatrix_index};
  pll_update_prob_matrices(partition, 0, matrix_indices, branch_lengths, 3);

  // make sure the parent clv is properly initialized TODO not needed for agressive range comp
  fill_without(partition, parent->clv_index, range, 1.0);

  // update the partial
  update_partial_ranged(partition, tree, range);

}

/* Adapted from LLPLL: pll_optimize.c
 * Original author: Diego Darriba
 */
static double recomp_iterative_ranged (pll_newton_tree_params_t * params,
                                double *best_lnl,
                                int radius,
                                Range& range)
{
  double lnl = 0.0,
     new_lnl = 0.0;
  pll_utree_t *tr_p, *tr_q, *tr_z;

  lnl = *best_lnl;
  tr_p = params->tree;
  tr_q = params->tree->next;
  tr_z = tr_q?tr_q->next:NULL;

  /* set Branch Length */
  assert(d_equals(tr_p->length, tr_p->back->length));

  double xmin, xguess, xmax;

  pll_update_sumtable (params->partition,
                       tr_p->clv_index,
                       tr_p->back->clv_index,
                       params->params_index,
                       params->freqs_index,
                       params->sumtable);

  xmin   = PLL_OPT_MIN_BRANCH_LEN + PLL_LBFGSB_ERROR;
  xmax   = PLL_OPT_MAX_BRANCH_LEN;
  xguess = tr_p->length;
  if (xguess < xmin || xguess > xmax)
    xguess = PLL_OPT_DEFAULT_BRANCH_LEN;

  double xres = pll_minimize_newton (xmin, xguess, xmax,
                                     10, &new_lnl, params,
                                     utree_derivative_func);

  if (pll_errno)
  {
    return PLL_FAILURE;
  }

  if (new_lnl >= 0)
    new_lnl = *best_lnl;

  /* ensure that new_lnl is not NaN */
  assert (new_lnl == new_lnl);

  if (new_lnl > lnl)
  {
    /* consolidate */
    tr_p->length = xres;
    tr_p->back->length = tr_p->length;

    pll_update_prob_matrices(params->partition,
                             params->params_index,
                             &(params->tree->pmatrix_index),
                             &(tr_p->length),1);

    *best_lnl = new_lnl;
  }
  else
  {
    /* revert */
   pll_update_prob_matrices(params->partition,
                             params->params_index,
                             &tr_p->pmatrix_index,
                             &tr_p->length, 1);

  }

  /* update children */
  if (radius && tr_q && tr_z)
  {
    /* update children 'Q'
     * CLV at P is recomputed with children P->back and Z->back
     * Scaler is updated by subtracting Q->back and adding P->back
     */

    update_partials_and_scalers(params->partition,
                      tr_q,
                      tr_p,
                      tr_z);

    /* eval */
    pll_newton_tree_params_t params_cpy;
    memcpy(&params_cpy, params, sizeof(pll_newton_tree_params_t));
    params_cpy.tree = tr_q->back;
    lnl = recomp_iterative_ranged (&params_cpy, best_lnl, radius-1, range);

    /* update children 'Z'
     * CLV at P is recomputed with children P->back and Q->back
     * Scaler is updated by subtracting Z->back and adding Q->back
     */

    update_partials_and_scalers(params->partition,
                      tr_z,
                      tr_q,
                      tr_p);

   /* eval */
    params_cpy.tree = tr_z->back;
    lnl = recomp_iterative_ranged (&params_cpy, best_lnl, radius-1, range);

    /* reset to initial state
     * CLV at P is recomputed with children Q->back and Z->back
     * Scaler is updated by subtracting P->back and adding Z->back
     */

    update_partials_and_scalers(params->partition,
                      tr_p,
                      tr_z,
                      tr_q);
  }

  return lnl;
}

/* Adapted from LLPLL: pll_optimize.c
 * Original author: Diego Darriba
 */
static double optimize_branch_triplet_newton(pll_partition_t * partition, pll_utree_t * tree, Range& range,
  const double tolerance, const int smoothings)
{
  unsigned int iters;
  double lnl = 0.0;

  /* get the initial likelihood score */
  lnl = compute_edge_logl_ranged(partition, tree, range);

  /**
   * preconditions:
   *    (1) CLVs must be updated towards 'tree'
   *    (2) Pmatrix indices must be **unique** for each branch
   */

  /* set parameters for N-R optimization */
  pll_newton_tree_params_t params;
  params.partition    = partition;
  params.tree         = tree;
  params.params_index = 0; // TODO change if we ever have multiple freqs or params
  params.freqs_index  = 0;
  params.sumtable     = 0;

  /* allocate the sumtable */
  if ((params.sumtable = (double *) pll_aligned_alloc(
      partition->sites * partition->rate_cats * partition->states *
      sizeof(double), partition->alignment)) == NULL)
  {
    throw runtime_error{"Cannot allocate memory for bl opt variables"};
  }

  iters = smoothings;
  while (iters)
  {
    double new_lnl = lnl;

    /* iterate on first edge */
    recomp_iterative (&params, &new_lnl, radius, keep_update);
    assert(new_lnl >= lnl);

    /* iterate on second edge */
    params.tree = params.tree->back;
    recomp_iterative (&params, &new_lnl, radius-1, keep_update);
    assert(new_lnl >= lnl);

    lnl = new_lnl;
    iters --;

    /* check convergence */
    if (fabs (new_lnl - lnl) < tolerance) iters = 0;
  }

  free(params.sumtable);

  return lnl;

}

double optimize_branch_triplet(pll_partition_t * partition, pll_utree_t * tree, Range& range)
{
  if (!tree->next)
    tree = tree->back;

  update_matrices_partial_triplet(partition, tree, range);

  auto logl = compute_edge_logl_ranged(partition, tree, range);

  auto cur_logl = -numeric_limits<double>::infinity();
  int smoothings = 8;

  while (fabs (cur_logl - logl) > OPT_EPSILON)
  {
    logl = cur_logl;
    cur_logl = optimize_branch_triplet_newton (partition,
                                              tree,
                                              OPT_BRANCH_EPSILON,
                                              smoothings);
    // re-update after change
    update_matrices_partial_triplet(partition, tree, range);
  }

  cur_logl = compute_edge_logl_ranged(partition, tree, range);

  return cur_logl;
}

static void traverse_update_partials(pll_utree_t * tree, pll_partition_t * partition,
    pll_utree_t ** travbuffer, double * branch_lengths, unsigned int * matrix_indices,
    pll_operation_t * operations)
{
  unsigned int num_matrices, num_ops;
  /* perform a full traversal*/
  assert(tree->next != nullptr);
  unsigned int traversal_size;
  // TODO this only needs to be done once, outside of this func. pass traversal size also
  // however this is practically nonexistent impact compared to clv comp
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

static double optimize_branch_lengths(pll_utree_t * tree, pll_partition_t * partition, pll_optimize_options_t& params,
  pll_utree_t ** travbuffer, double cur_logl, double lnl_monitor, int* smoothings)
{
  if (!tree->next)
    tree = tree->back;

  traverse_update_partials(tree, partition, travbuffer, params.lk_params.branch_lengths,
    params.lk_params.matrix_indices, params.lk_params.operations);

  pll_errno = 0; // hotfix

  cur_logl = -1 * pll_optimize_branch_lengths_iterative(
      partition, tree, 0, 0,
      OPT_BRANCH_EPSILON, *smoothings, 1); // last param = 1 means branch lengths are iteratively
                                      // updated during the call
  if (cur_logl+1e-6 < lnl_monitor)
    throw runtime_error{string("cur_logl < lnl_monitor: ") + to_string(cur_logl) + string(" : ")
    + to_string(lnl_monitor)};

  // reupdate the indices as they may have changed
  params.lk_params.where.unrooted_t.parent_clv_index = tree->clv_index;
  params.lk_params.where.unrooted_t.parent_scaler_index = tree->scaler_index;
  params.lk_params.where.unrooted_t.child_clv_index = tree->back->clv_index;
  params.lk_params.where.unrooted_t.child_scaler_index = tree->back->scaler_index;
  params.lk_params.where.unrooted_t.edge_pmatrix_index = tree->pmatrix_index;

  traverse_update_partials(tree, partition, travbuffer, params.lk_params.branch_lengths,
    params.lk_params.matrix_indices, params.lk_params.operations);
    cur_logl = pll_compute_edge_loglikelihood (partition, tree->clv_index,
                                                  tree->scaler_index,
                                                  tree->back->clv_index,
                                                  tree->back->scaler_index,
                                                  tree->pmatrix_index, 0);

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
  auto cur_logl = pll_compute_edge_loglikelihood (partition, tree->clv_index,
                                                tree->scaler_index,
                                                tree->back->clv_index,
                                                tree->back->scaler_index,
                                                tree->pmatrix_index, 0);

  // double cur_logl = -numeric_limits<double>::infinity();
  int smoothings = 32;
  double lnl_monitor = cur_logl;

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
  params.factr = OPT_FACTR;
  params.pgtol = OPT_PARAM_EPSILON;

  vector<pll_utree_t*> branches(nums.branches);
  auto num_traversed = utree_query_branches(tree, &branches[0]);
  assert (num_traversed == nums.branches);
  unsigned int branch_index = 0;

  double logl = cur_logl;

  if (opt_branches)
  {
    smoothings = 8;
    cur_logl = optimize_branch_lengths(branches[branch_index],
      partition, params, &travbuffer[0], cur_logl, lnl_monitor, &smoothings);

    lgr << "after hidden blo crouching tiger: " << to_string(cur_logl) << "\n";

  }

  #define RATE_MIN     1e-4
  #define RATE_MAX     1000000.
  double min_rates[6] = {RATE_MIN,RATE_MIN,RATE_MIN,RATE_MIN,RATE_MIN,RATE_MIN};
  double max_rates[6] = {RATE_MAX,RATE_MAX,RATE_MAX,RATE_MAX,RATE_MAX,RATE_MAX};


  do
  {
    branch_index = rand () % num_traversed;

    lgr << "Start: " << to_string(cur_logl) << "\n";
    logl = cur_logl;

    if (opt_model)
    {

      params.which_parameters = PLL_PARAMETER_SUBST_RATES;
      cur_logl = -pll_optimize_parameters_multidim(&params, min_rates, max_rates);

      lgr << "after rates: " << to_string(cur_logl) << "\n";

      if (opt_branches)
      {
        smoothings = 3;
        cur_logl = optimize_branch_lengths(branches[branch_index],
          partition, params, &travbuffer[0], cur_logl, lnl_monitor, &smoothings);

        lgr << "after blo 1: " << to_string(cur_logl) << "\n";

      }

      // params.which_parameters = PLL_PARAMETER_FREQUENCIES;
      // pll_optimize_parameters_multidim(&params, nullptr, nullptr);

      if (opt_branches)
      {
        smoothings = 3;
        cur_logl = optimize_branch_lengths(branches[branch_index],
          partition, params, &travbuffer[0], cur_logl, lnl_monitor, &smoothings);

        lgr << "after blo 2: " << to_string(cur_logl) << "\n";

      }

      // params.which_parameters = PLL_PARAMETER_PINV;
      // cur_logl = -1 * pll_optimize_parameters_brent(&params);
      params.which_parameters = PLL_PARAMETER_ALPHA;
      cur_logl = -pll_optimize_parameters_onedim(&params, 0.02, 10000.);

      lgr << "after alpha: " << to_string(cur_logl) << "\n";

    }

    if (opt_branches)
    {
      smoothings = 4;
      cur_logl = optimize_branch_lengths(branches[branch_index],
        partition, params, &travbuffer[0], cur_logl, lnl_monitor, &smoothings);

      lgr << "after blo 3: " << to_string(cur_logl) << "\n";

    }
  } while (fabs (cur_logl - logl) > OPT_EPSILON);

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
