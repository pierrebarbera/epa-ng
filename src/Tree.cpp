#include "Tree.h"

Tree::Tree(std::string* tree_file, std::string* msa_file)
{ 
  //detect filetypes
  //parse, build tree
  build_from_file(tree_file, &pll_utree_parse_newick);
  msa_from_file(msa_file)
}

Tree::~Tree()
{
	pll_partition_destroy(Tree::partition);
}

void Tree::build_tree_from_file(std::string* tree_file, pll_utree_t * (*tree_parse_f) (const char*, int*))
{
  int num_tip_nodes, num_nodes, num_branches, num_inner_nodes;

  /* first we call the appropriate pll parsing function to obtain a pll_utree structure, 
    on which our partition object will be based */
  auto tree = tree_parse_f(*tree_file, &num_tip_nodes);
  assert(num_tip_nodes == Tree::num_tip_nodes);

  set_missing_branch_length(tree, DEFAULT_BRANCH_LENGTH);

  // we then derive some numbers about the graph
  num_inner_nodes = num_tip_nodes - 2;
  num_nodes = num_inner_nodes + num_tip_nodes;
  num_branches = num_nodes - 1;

  // next, we obtain pointers to all tip nodes
  auto tip_nodes = (pll_utree_t  **)calloc(num_tip_nodes, sizeof(pll_utree_t *));
  pll_utree_query_tipnodes(tree, tip_nodes);

  /* allocate space to store sequences and headers */
  char ** headers = (char **)calloc(num_tip_nodes, sizeof(char *));
  char ** seqdata = (char **)calloc(num_tip_nodes, sizeof(char *));




  Tree::partition = pll_partition_create(num_tip_nodes,
                                   num_inner_nodes,
                                   STATES,
                                   sites,
                                   1,
                                   num_branches,
                                   RATE_CATS,
                                   num_inner_nodes,
                                   PLL_ATTRIB_ARCH_SSE);
}

// assumes msa file contains the correct number of reference sequences
// assumes references come first, then query
// assumes counting from 0
/* reads sequences from msa_file, returns them in <TODO usable format>
  msa_file: string specifying the file path
  stride: number of sequences to be read
*/
void* read_msa(std::string* msa_file, int stride)
{

  /* open the file */
  pll_fasta_t * file = pll_fasta_open(msa_file, pll_map_fasta);
  if (!file)
    error("Cannot open file %s", msa_file);

  char * sequence = NULL;
  char * header = NULL;
  long sequence_length;
  long header_length;
  long sequence_number;

  // skip ahead to

  /* read sequences and make sure they are all of the same length */
  int num_sites = NULL;
  for (i = first; pll_fasta_getnext(fp, &header, &header_length, &sequence,
                                &sequence_length, &sequence_number) && (i < last); ++i)
  {

    if (!num_sites && num_sites != sequence_length)
      error("MSA file does not contain equal size sequences\n");

    if (!num_sites) num_sites = sequence_length;

    headers[i] = header;
    seqdata[i] = sequence;
  }

  /* did we stop reading the file because we reached EOF? */
  if (pll_errno != PLL_ERROR_FILE_EOF)
    error("Error while reading file " +  msa_file);

  pll_fasta_close(file);

  if (num_sites == -1)
    error("Unable to read alignment");

  if (i != tip_nodes_count)
    error("Some taxa are missing from MSA file");

}

void* read_msa(void* fp, int stride) 
{

}
