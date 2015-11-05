#include "Tree.h"

Tree::Tree(const std::string& tree_file, const std::string& msa_file)
{ 
  //TODO handle filetypes other than newick/fasta

  //parse, build tree
  Tree::msa = new MSA(msa_file);
  build_partition_from_file(tree_file, &pll_utree_parse_newick);
}

Tree::~Tree()
{
	pll_partition_destroy(Tree::partition);
  delete this->msa;
}

void Tree::build_partition_from_file(const std::string& tree_file, pll_utree_t * (*tree_parse_f) (const char*, int*))
{
  int num_tip_nodes, num_nodes, num_branches, num_inner_nodes;

  /* first we call the appropriate pll parsing function to obtain a pll_utree structure, 
    on which our partition object will be based */
  auto tree = tree_parse_f(tree_file.c_str(), &num_tip_nodes);
  assert(num_tip_nodes == Tree::num_tip_nodes);

  set_missing_branch_length(tree, DEFAULT_BRANCH_LENGTH);

  // we then derive some numbers about the graph
  num_inner_nodes = num_tip_nodes - 2;
  num_nodes = num_inner_nodes + num_tip_nodes;
  num_branches = num_nodes - 1;

  // next, we obtain pointers to all tip nodes
  auto tip_nodes = (pll_utree_t  **)calloc(num_tip_nodes, sizeof(pll_utree_t *));
  pll_utree_query_tipnodes(tree, tip_nodes);


  Tree::partition = pll_partition_create(num_tip_nodes,
                                   num_inner_nodes,
                                   STATES,
                                   Tree::msa->num_sites,
                                   1,
                                   num_branches,
                                   RATE_CATS,
                                   num_inner_nodes,
                                   PLL_ATTRIB_ARCH_SSE);

  /* initialize the array of base frequencies */
  double frequencies[4] = { 0.17, 0.19, 0.25, 0.39 };

  /* substitution rates for the 4x4 GTR model. This means we need exactly
     (4*4-4)/2 = 6 values, i.e. the number of elements above the diagonal */
  double subst_params[6] = {1,1,1,1,1,1};

  /* we'll use 4 rate categories, and currently initialize them to 0 */
  double rate_cats[4] = {0};

  /* compute the discretized category rates from a gamma distribution
     with alpha shape 1 and store them in rate_cats  */
  pll_compute_gamma_cats(1, 4, rate_cats);

  /* set frequencies at model with index 0 (we currently have only one model) */
  pll_set_frequencies(partition, 0, frequencies);

  /* set 6 substitution parameters at model with index 0 */
  pll_set_subst_params(partition, 0, subst_params);

  /* set rate categories */
  pll_set_category_rates(partition, rate_cats);

}


