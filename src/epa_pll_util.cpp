#include "epa_pll_util.hpp"

#include <unordered_map>
#include <string>
#include <algorithm>
#include <cmath>
#include <cstdlib>

#include "pll_util.hpp"

#define OPT_EPSILON       1e-1
#define OPT_PARAM_EPSILON 1e-1

/* if set, the parameters are no longer optimized when
 * they do not improve the likelihood at one iteration */
#define CHECK_LOCAL_CONVERGENCE 1

using namespace std;

void link_tree_msa(pll_utree_t * tree, pll_partition_t * partition,
              const MSA& msa, const unsigned int num_tip_nodes)
{
  // obtain pointers to all tip nodes
  vector<pll_utree_t*> tip_nodes(num_tip_nodes);
  pll_utree_query_tipnodes(tree, &tip_nodes[0]);

  // and associate the sequences from the MSA file with the correct tips
  /* create a hash table of size num_tip_nodes */
  unordered_map<string, unsigned int> map; // mapping labels to tip clv indices

  /* populate the hash table with tree tip labels */
  for (unsigned int i = 0; i < num_tip_nodes; ++i)
    map[tip_nodes[i]->label] = i;

  /* find sequences in hash table and link them with the corresponding taxa */
  for (auto const &s : msa)
  {
    auto map_value = map.find(s.header());

    if (map_value == map.end())
      throw runtime_error{string("Sequence with header does not appear in the tree: ") + s.header()};

    // associates the sequence with the tip by calculating the tips clv buffers
    pll_set_tip_states(partition, map_value->second, pll_map_nt, s.sequence().c_str());
  }
};

void precompute_clvs(pll_utree_t * tree, pll_partition_t * partition, const Tree_Numbers& nums)
{
  unsigned int num_matrices, num_ops;

  /* various buffers for creating a postorder traversal and operations structures */
  vector<pll_utree_t*> travbuffer(nums.nodes);
  vector<double> branch_lengths(nums.branches);
  vector<unsigned int> matrix_indices(nums.branches);
  vector<pll_operation_t> operations(nums.nodes);

  // get a list of all tip nodes
  vector<pll_utree_t*> tip_nodes(nums.tip_nodes);
  pll_utree_query_tipnodes(tree, &tip_nodes[0]);

  /* adjust clv indices such that every direction has its own */
  set_unique_clv_indices(tree, nums.tip_nodes);

  for (auto node : tip_nodes)
  {
    /* perform a partial postorder traversal of the unrooted tree  starting at the current tip
      and returning every node whose clv in the direction of the tip hasn't been calculated yet*/
    unsigned int traversal_size;
    if(pll_utree_traverse(node->back, cb_partial_traversal, &travbuffer[0], &traversal_size)
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
  }
}

void split_combined_msa(MSA& source, MSA& target, pll_utree_t * tree, unsigned int num_tip_nodes)
{
  vector<pll_utree_t*> tip_nodes(num_tip_nodes);
  pll_utree_query_tipnodes(tree, &tip_nodes[0]);

  auto falsegroup_begin = partition(source.begin(), source.end(),
    [&tip_nodes](const Sequence& em)
    {
      return find(tip_nodes.begin(), tip_nodes.end(), em) != tip_nodes.end();
    });
  target.num_sites(source.num_sites());
  target.move_sequences(falsegroup_begin, source.end());
  source.erase(falsegroup_begin, source.end());
}

bool operator==(const pll_utree_t * node, const Sequence& s)
{
  return s.header().compare(node->label) == 0;
}

bool operator==(const Sequence& s, const pll_utree_t * node)
{
  return operator==(node, s);
}

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
