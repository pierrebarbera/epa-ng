#include "optimize.hpp"

#include <vector>
#include <cmath>
#include <cstdlib>
#include <stdexcept>

#include "pll_util.hpp"

using namespace std;

void optimize(pll_utree_t * tree, pll_partition_t * partition, const Tree_Numbers &nums,
    Model& model)
{
  unsigned int num_matrices, num_ops;

  /* various buffers for creating a postorder traversal and operations structures */
  vector<pll_utree_t*> travbuffer(nums.nodes);
  vector<double> branch_lengths(nums.branches);
  vector<unsigned int> matrix_indices(nums.branches);
  vector<pll_operation_t> operations(nums.nodes);

  vector<pll_utree_t *> inner_nodes(nums.inner_nodes);
  pll_utree_query_innernodes(tree, &inner_nodes[0]);

  /* perform a full traversal*/
  unsigned int traversal_size;
  if(pll_utree_traverse(tree, cb_full_traversal, &travbuffer[0], &traversal_size)
              != PLL_SUCCESS)
    throw runtime_error{"Function pll_utree_traverse() requires inner nodes as parameters"};

  /* given the computed traversal descriptor, generate the operations
     structure, and the corresponding probability matrix indices that
     may need recomputing */
  pll_utree_create_operations(&travbuffer[0],
                              traversal_size,
                              &branch_lengths[0],
                              &matrix_indices[0],
                              &operations[0],
                              &num_matrices,
                              &num_ops);

  pll_update_prob_matrices(partition,
                           0,             // use model 0
                           &matrix_indices[0],// matrices to update
                           &branch_lengths[0],
                           num_matrices); // how many should be updated

  /* use the operations array to compute all num_ops inner CLVs. Operations
     will be carried out sequentially starting from operation 0 towrds num_ops-1 */
  pll_update_partials(partition, &operations[0], num_ops);

  // compute logl once to give us a logl starting point
  auto logl = pll_compute_edge_loglikelihood (partition, tree->clv_index,
                                                tree->scaler_index,
                                                tree->back->clv_index,
                                                tree->back->scaler_index,
                                                tree->pmatrix_index, 0);

  pll_optimize_options_t params;
  params.lk_params.partition = partition;
  params.lk_params.operations = &operations[0];
  params.lk_params.branch_lengths = &branch_lengths[0];
  params.lk_params.matrix_indices = &matrix_indices[0];
  params.lk_params.alpha_value = 0;
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
  params.subst_params_symmetries = &(model.symmetries())[0];
  params.factr = 1e7;
  params.pgtol = OPT_PARAM_EPSILON;

  params.which_parameters = PLL_PARAMETER_BRANCHES_SINGLE;

  double cur_logl = logl - 10;
  int smoothings = 1;
  double lnl_monitor = logl;

  while (fabs (cur_logl - logl) > OPT_EPSILON)
  {
    logl = cur_logl;

    // set current node to a random inner node to avoid local maxima
    int inner_index = rand () % nums.inner_nodes;
    tree = inner_nodes[inner_index];

    pll_utree_traverse(tree, cb_full_traversal, &travbuffer[0], &traversal_size);

    /* given the computed traversal descriptor, generate the operations
       structure, and the corresponding probability matrix indices that
       may need recomputing */
    pll_utree_create_operations(&travbuffer[0],
                                traversal_size,
                                &branch_lengths[0],
                                &matrix_indices[0],
                                &operations[0],
                                &num_matrices,
                                &num_ops);

    pll_update_prob_matrices(partition,
                             0,             // use model 0
                             &matrix_indices[0],// matrices to update
                             &branch_lengths[0],
                             num_matrices); // how many should be updated

    pll_update_partials(partition, &operations[0], num_ops);

    // reupdate params for current node
    params.lk_params.where.unrooted_t.parent_clv_index = tree->clv_index;
    params.lk_params.where.unrooted_t.parent_scaler_index =
        tree->scaler_index;
    params.lk_params.where.unrooted_t.child_clv_index = tree->back->clv_index;
    params.lk_params.where.unrooted_t.child_scaler_index =
        tree->back->scaler_index;
    params.lk_params.where.unrooted_t.edge_pmatrix_index =
        tree->pmatrix_index;

    pll_errno = 0; // hotfix

    cur_logl = -1 * pll_optimize_branch_lengths_iterative (
        partition, tree, params.params_index, params.lk_params.freqs_index,
        params.pgtol, smoothings++, 1); // last param = 1 means branch lengths are iteratively
                                        // updated during the call

    assert(cur_logl >= lnl_monitor);

    if (CHECK_LOCAL_CONVERGENCE && fabs(cur_logl - lnl_monitor) < OPT_EPSILON)
      break;

    lnl_monitor = cur_logl;

  }
}
