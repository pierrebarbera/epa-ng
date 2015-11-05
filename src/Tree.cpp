#include "Tree.h"

using namespace std;

Tree::Tree(const string& tree_file, const string& msa_file, Model& model) : model(model)
{ 

  //TODO handle filetypes other than newick/fasta

  //parse, build tree
  msa = new MSA(msa_file);
  build_partition_from_file(tree_file);
}

Tree::~Tree()
{
	pll_partition_destroy(Tree::partition);
  delete msa;
}


//TODO too long, refactor
void Tree::build_partition_from_file(const string& tree_file)
{
  int num_tip_nodes, num_nodes, num_branches, num_inner_nodes;

  /* first we call the appropriate pll parsing function to obtain a pll_utree structure, 
    on which our partition object will be based */
  auto tree = pll_utree_parse_newick(tree_file.c_str(), &num_tip_nodes);

  if (num_tip_nodes < 3)
    throw runtime_error{"Number of tip nodes too small"};

  set_missing_branch_length(tree, DEFAULT_BRANCH_LENGTH);

  // we then derive some numbers about the graph
  num_inner_nodes = num_tip_nodes - 2;
  num_nodes = num_inner_nodes + num_tip_nodes;
  num_branches = num_nodes - 1;

  // next, we obtain pointers to all tip nodes
  auto tip_nodes = new pll_utree_t*[num_tip_nodes];
  pll_utree_query_tipnodes(tree, tip_nodes);

  /* create a libc hash table of size num_tip_nodes */
  int tmp = num_tip_nodes;
  hcreate(tmp);

  /* populate a libc hash table with tree tip labels */
  auto data = new int[num_tip_nodes];
  for (int i = 0; i < num_tip_nodes; ++i)
  {
    data[i] = i;
    ENTRY entry;
    entry.key = tip_nodes[i]->label;
    entry.data = (void *)(data+i);
    hsearch(entry, ENTER);
  }


  Tree::partition = pll_partition_create(num_tip_nodes,
                                   num_inner_nodes,
                                   STATES,
                                   msa->num_sites,
                                   1,
                                   num_branches,
                                   RATE_CATS,
                                   num_inner_nodes,
                                   PLL_ATTRIB_ARCH_SSE);

  double rate_cats[RATE_CATS] = {0};

  /* compute the discretized category rates from a gamma distribution
     with alpha shape 1 and store them in rate_cats  */
  // TODO gamma dist is input?
  pll_compute_gamma_cats(1, 4, rate_cats);

  /* set frequencies at model with index 0 */
  pll_set_frequencies(partition, 0, &(model.base_frequencies[0]));

  /* set substitution parameters at model with index 0 */
  pll_set_subst_params(partition, 0, &(model.substitution_rates[0]));

  /* set rate categories */
  pll_set_category_rates(partition, rate_cats);

  /* find sequences in hash table and link them with the corresponding taxa */
  for (int i = 0; i < num_tip_nodes; ++i)
  {
    string header, sequence;
    tie(header, sequence) = msa->get(i);
    ENTRY query;
    query.key = &header[0];
    ENTRY * found = NULL;

    found = hsearch(query,FIND);

    if (!found)
      throw runtime_error{string("Sequence with header does not appear in the tree: ") + header};
        
    int tip_clv_index = *((int *)(found->data));

    pll_set_tip_states(partition, tip_clv_index, pll_map_nt, sequence.c_str());
  }

  /* destroy hash table and free related data */
  hdestroy();
  delete [] tip_nodes;
  delete [] data;

  //TODO reference part of the MSA can now be freed
}

void precompute_clvs()
{

}

