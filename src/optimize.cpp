#include "optimize.hpp"

#include <vector>
#include <cmath>
#include <cstdlib>
#include <cassert>
#include <stdexcept>

#include "pll_util.hpp"
#include "constants.hpp"

using namespace std;

void optimize(pll_utree_t * tree, pll_partition_t * partition, const Tree_Numbers &nums)
{
  unsigned int num_matrices, num_ops;

  /* various buffers for creating a postorder traversal and operations structures */
  vector<pll_utree_t*> travbuffer(nums.nodes);
  vector<double> branch_lengths(nums.branches);
  vector<unsigned int> matrix_indices(nums.branches);
  vector<pll_operation_t> operations(nums.nodes);

  vector<pll_utree_t *> branches(nums.branches);
  utree_query_branches(tree, &branches[0]);

  /* perform a full traversal*/
  assert(tree->next != nullptr);
  unsigned int traversal_size;
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

  /* use the operations array to compute all num_ops inner CLVs. Operations
     will be carried out sequentially starting from operation 0 towrds num_ops-1 */
  pll_update_partials(partition, &operations[0], num_ops);

  // compute logl once to give us a logl starting point
  auto logl = pll_compute_edge_loglikelihood (partition, tree->clv_index,
                                                tree->scaler_index,
                                                tree->back->clv_index,
                                                tree->back->scaler_index,
                                                tree->pmatrix_index, 0);

  double cur_logl = -99999999999999.0;
  int smoothings = 1;
  double lnl_monitor = logl;

  while (fabs (cur_logl - logl) > OPT_EPSILON)
  {
    logl = cur_logl;

    // set current node to a random inner node to avoid local maxima
    tree = branches[rand () % nums.branches];

    if (!tree->next)
      tree = tree->back;

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
