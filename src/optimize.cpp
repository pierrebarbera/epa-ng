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

void optimize_branch_lengths(pll_utree_t * tree, pll_partition_t * partition, const Tree_Numbers &nums)
{
  /* various buffers for creating a postorder traversal and operations structures */
  vector<pll_utree_t*> travbuffer(nums.nodes);
  vector<double> branch_lengths(nums.branches);
  vector<unsigned int> matrix_indices(nums.branches);
  vector<pll_operation_t> operations(nums.nodes);

  vector<pll_utree_t *> branches(nums.branches);
  utree_query_branches(tree, &branches[0]);

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

  while (fabs (cur_logl - logl) > OPT_EPSILON)
  {
    logl = cur_logl;

    // set current node to a random inner node to avoid local maxima
    tree = branches[rand () % nums.branches];

    if (!tree->next)
      tree = tree->back;

    traverse_update_partials(tree, partition, &travbuffer[0], &branch_lengths[0],
      &matrix_indices[0], &operations[0]);

    pll_errno = 0; // hotfix

    cur_logl = -1 * pll_optimize_branch_lengths_iterative (
        partition, tree, 0, 0,
        OPT_PARAM_EPSILON, smoothings++, 1); // last param = 1 means branch lengths are iteratively
                                        // updated during the call
    if(cur_logl+1e-6 < lnl_monitor)
      throw runtime_error{string("cur_logl < lnl_monitor: ") + to_string(cur_logl) + string(" : ")
      + to_string(lnl_monitor)};


    if (CHECK_LOCAL_CONVERGENCE && fabs(cur_logl - lnl_monitor) < OPT_EPSILON)
      break;

    lnl_monitor = cur_logl;

  }
}

void optimize_model_params(Model& model, pll_utree_t * tree, pll_partition_t * partition,
  const Tree_Numbers& nums)
{
  auto symmetries = (&(model.symmetries())[0]);

  // sadly we explicitly need these buffers here and in the params structure
  vector<pll_utree_t*> travbuffer(nums.nodes);
  vector<double> branch_lengths(nums.branches);
  vector<unsigned int> matrix_indices(nums.branches);
  vector<pll_operation_t> operations(nums.nodes);

  traverse_update_partials(tree, partition, &travbuffer[0], &branch_lengths[0],
    &matrix_indices[0], &operations[0]);

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

  // params.which_parameters = PLL_PARAMETER_ALPHA;
  // pll_optimize_parameters_brent(&params);
  //
  // params.which_parameters = PLL_PARAMETER_SUBST_RATES;
  // pll_optimize_parameters_lbfgsb(&params);
  //
  // params.which_parameters = PLL_PARAMETER_FREQUENCIES;
  // pll_optimize_parameters_lbfgsb(&params);
  //
  // params.which_parameters = PLL_PARAMETER_PINV;
  // pll_optimize_parameters_lbfgsb(&params);

  params.which_parameters = PLL_PARAMETER_FREQUENCIES
  | PLL_PARAMETER_ALPHA | PLL_PARAMETER_SUBST_RATES;// | PLL_PARAMETER_PINV;
  pll_optimize_parameters_lbfgsb(&params);

  // params.which_parameters = PLL_PARAMETER_BRANCHES_ITERATIVE;
  // pll_optimize_parameters_lbfgsb(&params);

  // update epa model object as well
  model.alpha(params.lk_params.alpha_value);
  model.substitution_rates(partition->subst_params[0], 6);
  model.base_frequencies(partition->frequencies[params.params_index], partition->states);
}

void compute_and_set_empirical_frequencies(pll_partition_t * partition)
{
  double * empirical_freqs = pll_compute_empirical_frequencies (partition);
  pll_set_frequencies (partition, 0, 0, empirical_freqs);
  free (empirical_freqs);
}
