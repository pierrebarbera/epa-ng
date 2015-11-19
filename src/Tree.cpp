#include "Tree.hpp"

#include <stdexcept>

#include "pll_util.hpp"
#include "epa_pll_util.hpp"
#include "file_io.hpp"
#include "Sequence.hpp"

using namespace std;

Tree::Tree(const string& tree_file, const string& msa_file, const Model& model) : model_(model)
{

  // TODO handle filetypes other than newick/fasta


  //parse, build tree
  ref_msa_ = build_MSA_from_file(msa_file);

  build_partition_from_file(tree_file);

  link_tree_msa(tree_, partition_, ref_msa_, num_tip_nodes_);

  precompute_clvs();
}

Tree::~Tree()
{
  // free data segment of tree nodes
  utree_free_node_data(tree_);

	pll_partition_destroy(partition_);

}

void Tree::build_partition_from_file(const string& tree_file)
{
  int  num_nodes, num_inner_nodes;

  /* first we call the appropriate pll parsing function to obtain a pll_utree structure,
    on which our partition object will be based */
  tree_ = pll_utree_parse_newick(tree_file.c_str(), &num_tip_nodes_);

  if (num_tip_nodes_ < 3)
    throw runtime_error{"Number of tip nodes too small"};

  set_missing_branch_length(tree_, DEFAULT_BRANCH_LENGTH);

  // we then derive some numbers about the graph
  num_inner_nodes = num_tip_nodes_ - 2;
  num_nodes = num_inner_nodes + num_tip_nodes_;
  num_branches_ = num_nodes - 1;

  partition_ = pll_partition_create(num_tip_nodes_,
                                   num_inner_nodes * 3, //number of extra clv buffers: 3 for every direction on the node
                                   STATES,
                                   ref_msa_.num_sites(),
                                   1,
                                   num_branches_,
                                   RATE_CATS,
                                   num_inner_nodes * 3, /* number of extra scaler buffers */
                                   PLL_ATTRIB_ARCH_SSE);

  double rate_cats[RATE_CATS] = {0};

  /* compute the discretized category rates from a gamma distribution
     with alpha shape 1 and store them in rate_cats  */
  // TODO gamma dist is input?
  pll_compute_gamma_cats(1, 4, rate_cats);

  /* set frequencies at model with index 0 */
  pll_set_frequencies(partition_, 0, &(model_.base_frequencies()[0]));

  /* set substitution parameters at model with index 0 */
  pll_set_subst_params(partition_, 0, &(model_.substitution_rates()[0]));

  /* set rate categories */
  pll_set_category_rates(partition_, rate_cats);

  // TODO reference part of the MSA can now be freed
}

void Tree::precompute_clvs()
{
  //rederive the numbers
  auto num_inner_nodes = num_tip_nodes_ - 2;
  auto num_nodes = num_inner_nodes + num_tip_nodes_;
  //num_branches_ = num_nodes - 1;

  int num_matrices, num_ops;

  /* buffer for creating a postorder traversal structure */
  auto travbuffer = new pll_utree_t*[num_nodes];
  auto branch_lengths = new double[num_branches_];
  auto matrix_indices = new int[num_branches_];
  auto operations = new pll_operation_t[num_nodes];
  // TODO num_nodes too many? whats the upper bound considering tip traversal?

  // get a list of all tip nodes
  // TODO refactor candidate
  auto tip_nodes = new pll_utree_t*[num_tip_nodes_];
  pll_utree_query_tipnodes(tree_, tip_nodes);

  /* adjust clv indices such that every direction has its own */
  set_unique_clv_indices(tree_, num_tip_nodes_);

  for (int i = 0; i < num_tip_nodes_; ++i)
  {
    /* perform a partial postorder traversal of the unrooted tree */

    int traversal_size = pll_utree_traverse(tip_nodes[i]->back,
                                            cb_partial_traversal,
                                            travbuffer);
    if (traversal_size == -1)
      throw runtime_error{"Function pll_utree_traverse() requires inner nodes as parameters"};

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

    pll_update_prob_matrices(partition_,
                             0,             // use model 0
                             matrix_indices,// matrices to update
                             branch_lengths,
                             num_matrices); // how many should be updated

    /* use the operations array to compute all num_ops inner CLVs. Operations
       will be carried out sequentially starting from operation 0 towrds num_ops-1 */
    pll_update_partials(partition_, operations, num_ops);
  }

  delete [] travbuffer;
  delete [] tip_nodes;
  delete [] branch_lengths;
  delete [] matrix_indices;
  delete [] operations;

}

void Tree::place(const MSA &msa) const
{
  // get all edges
  auto node_list = new pll_utree_t*[num_branches_];
  auto num_traversed = utree_query_branches(tree_, node_list);

  // TODO create output class

  // place all s on every edge
  for (auto s : msa)
    for (int i = 0; i < num_traversed; ++i) {
      place_on_edge(s, node_list[i]);
    }
}

/* Compute loglikelihood of a Sequence, when placed on the edge defined by node */
double Tree::place_on_edge(Sequence& s, pll_utree_t * node) const
{
  int old_left_clv_index = 0;
  int old_right_clv_index = 1;
  int new_tip_clv_index = 2;
  int inner_clv_index = 3;
  pll_utree_t * old_left = node;
  pll_utree_t * old_right = node->back;

  /* create new tiny tree including the given nodes and a new node for the sequence s

               [2]
             new_tip
                |
                |
              inner [3]
             /     \
            /       \
      old_left     old_right
        [0]           [1]

    numbers in brackets are the nodes clv indices
  */

  // stack (TODO) new partition with 3 tips, 1 inner node
  //partition_create_stackd(tiny_partition,
  pll_partition_t * tiny_partition =
  pll_partition_create(   3, // tips
                          1, // extra clv's
                          partition_->states,
                          partition_->sites,
                          partition_->rate_matrices,
                          3, // number of prob. matrices (one per unique branch length)
                          partition_->rate_cats,
                          4, // number of scale buffers (one per node)
                          partition_->attributes);

  // shallow copy model params
  tiny_partition->rates = partition_->rates;
  tiny_partition->subst_params = partition_->subst_params;
  tiny_partition->frequencies = partition_->frequencies;

  // shallow copy 2 existing nodes clvs
  tiny_partition->clv[old_left_clv_index] =
                          partition_->clv[old_left->clv_index];
  tiny_partition->clv[old_right_clv_index] =
                          partition_->clv[old_right->clv_index];

  // shallow copy scalers TODO is this right?
  tiny_partition->scale_buffer[old_left_clv_index] =
                          partition_->scale_buffer[old_left->scaler_index];
  tiny_partition->scale_buffer[old_right_clv_index] =
                          partition_->scale_buffer[old_right->scaler_index];

  // init the new tip with s.sequence(), branch length
  pll_set_tip_states(tiny_partition, new_tip_clv_index, pll_map_nt, s.sequence().c_str());

  // set up branch lengths
  //int num_unique_branch_lengths = 2;
  /* heuristic insertion as described in EPA paper from 2011 (Berger et al.):
    original branch, now split by "inner", or base, node of the inserted sequence,
    defines the new branch lengths between inner and old left/right respectively
    as old branch length / 2.
    The new branch leading from inner to the new tip is initialized with length 0.9,
    which is the default branch length in RAxML.
    */
  double branch_lengths[2] = { old_left->length / 2, 0.9};
  int matrix_indices[2] = { 0, 1 };

  // TODO Newton-Raphson

  // use branch lengths to compute the probability matrices
  pll_update_prob_matrices(tiny_partition, 0, matrix_indices,
                      branch_lengths, 2);

  /* Creating a single operation for the inner node computation.
    Here we specify which clvs and pmatrices to use, so we don't need to mess with
    what the internal tree structure points to */
  pll_operation_t ops[1];
  ops[0].parent_clv_index    = inner_clv_index;
  ops[0].child1_clv_index    = old_left_clv_index;
  ops[0].child2_clv_index    = old_right_clv_index;
  ops[0].child1_matrix_index = 0;// TODO depends on NR vs heuristic
  ops[0].child2_matrix_index = 0;// TODO depends on NR vs heuristic
  ops[0].parent_scaler_index = PLL_SCALE_BUFFER_NONE;// TODO this should be inner_clv_index once scale buffers are fixed
  ops[0].child1_scaler_index = old_left_clv_index;
  ops[0].child2_scaler_index = old_right_clv_index;

  // use update_partials to compute the clv poining toward the new tip
  pll_update_partials(tiny_partition, ops, 1);

  // compute the loglikelihood using inner node and new tip
  double likelihood =  pll_compute_edge_loglikelihood(tiny_partition,
                                        new_tip_clv_index,
                                        new_tip_clv_index,// scaler_index
                                        inner_clv_index,
                                        inner_clv_index,  // scaler_index
                                        1,// matrix index of branch TODO depends on NR
                                        0);// freq index
  //pll_partition_destroy(tiny_partition);
  return likelihood;
}

void Tree::visit(std::function<void(pll_partition_t *, pll_utree_t *)> f)
{
  f(partition_, tree_);
}
