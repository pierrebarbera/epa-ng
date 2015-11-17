#include "Tree.hpp"

#include <stdexcept>
#include <search.h>

#include "pll_util.hpp"
#include "file_io.hpp"
#include "Sequence.hpp"

using namespace std;

Tree::Tree(const string& tree_file, const string& msa_file, Model& model) : model_(model)
{

  // TODO handle filetypes other than newick/fasta


  //parse, build tree
  ref_msa_ = build_MSA_from_file(msa_file);

  build_partition_from_file(tree_file);

  link_tree_msa();

  precompute_clvs();
}

Tree::~Tree()
{
  // free data segment of tree nodes
  free_node_data(tree_);

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

void Tree::link_tree_msa()
{
  // next, we obtain pointers to all tip nodes
  auto tip_nodes = new pll_utree_t*[num_tip_nodes_];
  pll_utree_query_tipnodes(tree_, tip_nodes);

  /*
  // and associate the sequences from the MSA file with the correct tips
  // TODO this is n^2. potential performance hog
  bool found = false;
  for (int i = 0; i < num_tip_nodes_; ++i)
  {
    string header, sequence;
    tie(header, sequence) = msa->get(i);

    for (int j = 0; j < num_tip_nodes_ && !found; ++j)
    {
      if(header.compare(tip_nodes[j]->label) == 0)
      {
        found = true;
        pll_set_tip_states(partition, j, pll_map_nt, sequence.c_str());
      }
    }
    if (!found)
      throw runtime_error{string("Sequence with header does not appear in the tree: ") + header};
    found = false;
  }*/

  /* create a libc hash table of size num_tip_nodes_ */
  hcreate(num_tip_nodes_);

  /* populate a libc hash table with tree tip labels */
  auto data = new int[num_tip_nodes_];
  for (int i = 0; i < num_tip_nodes_; ++i)
  {
    data[i] = i;
    ENTRY entry;
    entry.key = tip_nodes[i]->label;
    entry.data = (void *)(&data[i]);
    hsearch(entry, ENTER);
  }

  /* find sequences in hash table and link them with the corresponding taxa */
  for (auto s : ref_msa_)
  {
    ENTRY query;
    query.key = &(s.header()[0]); //to get char* from std::string
    ENTRY * found = NULL;

    found = hsearch(query,FIND);

    if (!found)
      throw runtime_error{string("Sequence with header does not appear in the tree: ") + s.header()};

    int tip_clv_index = *((int *)(found->data));

    pll_set_tip_states(partition_, tip_clv_index, pll_map_nt, s.sequence().c_str());
  }

  /* destroy hash table and free related data */
  hdestroy();
  delete [] tip_nodes;
  delete [] data;
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
  pll_utree * inner;
  pll_utree * new_tip;
  pll_utree * old_left = node;
  pll_utree * old_right = node->back;

  /* create new tiny tree including the given nodes and a new node for the sequence s */

  // stack new partition with 3 tips, 1 inner node
  pll_partition_t * tiny_partition =
  partition_create_stackd(3, // tips
                          1, // extra clv's
                          partition_->states,
                          partition_->sites,
                          partition_->rate_matrices,
                          3, // number of prob. matrices (one per unique branch length)
                          partition_->rate_cats,
                          4, // number of scale buffers (one per node)
                          partition_->attributes);

  // shallow copy 2 existing nodes

  // init the new tip with s.sequence(), branch length
  s.sequence();


  // init the new inner node with proper branch length indices

  // create a single operation for the inner node computation

  // use update_partials to compute the clv poinint toward the new tip

  // compute the loglikelihood using inner node and new tip
  return pll_compute_edge_loglikelihood(tiny_partition,
                                        inner->clv_index,
                                        inner->scaler_index,
                                        new_tip->clv_index,
                                        new_tip->scaler_index,
                                        0,// matrix index of branch
                                        0);// freq index
}

void Tree::visit(std::function<void(pll_partition_t *, pll_utree_t *)> f)
{
  f(partition_, tree_);
}
