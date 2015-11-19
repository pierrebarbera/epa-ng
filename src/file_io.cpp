#include "file_io.hpp"

#include <stdexcept>

#include "pllhead.hpp"
#include "MSA.hpp"
#include "pll_util.hpp"

using namespace std;

/* reads in sequences from a file
  msa_file: string specifying the file path
*/
MSA build_MSA_from_file(const string& msa_file)
{
  /* open the file */
  auto file = pll_fasta_open(msa_file.c_str(), pll_map_fasta);
  if (!file)
    throw runtime_error{string("Cannot open file") + msa_file};

  char * sequence = NULL;
  char * header = NULL;
  long sequence_length;
  long header_length;
  long sequence_number;

  /* read sequences and make sure they are all of the same length */
  int sites = 0;

  /* read the first sequence seperately, so that the MSA object can be constructed */
  pll_fasta_getnext(file, &header, &header_length, &sequence, &sequence_length, &sequence_number);
  sites = sequence_length;

  if (sites == -1 || sites == 0)
    throw runtime_error{"Unable to read MSA file"};

  auto msa = MSA(sites);
  msa.append(header, sequence);

  /* read the rest */
  while (pll_fasta_getnext(file, &header, &header_length, &sequence, &sequence_length, &sequence_number))
  {

    if (sites && sites != sequence_length)
      throw runtime_error{"MSA file does not contain equal size sequences"};

    if (!sites) sites = sequence_length;

    msa.append(header, sequence);
  }

  if (pll_errno != PLL_ERROR_FILE_EOF)
    throw runtime_error{string("Error while reading file:") +  msa_file};

  pll_fasta_close(file);

  return msa;
}

tuple<pll_partition_t *, pll_utree_t *>  build_partition_from_file(const string& tree_file, const Model& model,
                  Tree_Numbers& nums,  const int num_sites)
{
  unsigned int num_tip_nodes;

  /* first we call the appropriate pll parsing function to obtain a pll_utree structure,
    on which our partition object will be based */
  auto tree = pll_utree_parse_newick(tree_file.c_str(), &num_tip_nodes);

  // we then derive some numbers about the graph
  nums.init(num_tip_nodes);

  if (nums.tip_nodes < 3)
    throw runtime_error{"Number of tip nodes too small"};

  set_missing_branch_length(tree, DEFAULT_BRANCH_LENGTH);


  auto partition = pll_partition_create(nums.tip_nodes,
                                   nums.inner_nodes * 3, //number of extra clv buffers: 3 for every direction on the node
                                   STATES,
                                   num_sites,
                                   0, // number of mixture models
                                   1,
                                   nums.branches,
                                   RATE_CATS,
                                   nums.inner_nodes * 3, /* number of extra scaler buffers */
                                   PLL_ATTRIB_ARCH_SSE);

  double rate_cats[RATE_CATS] = {0};

  /* compute the discretized category rates from a gamma distribution
     with alpha shape 1 and store them in rate_cats  */
  // TODO gamma dist is input?
  pll_compute_gamma_cats(1, 4, rate_cats);

  /* set frequencies at model with index 0 */
  pll_set_frequencies(partition, 0, 0, &(model.base_frequencies()[0]));

  /* set substitution parameters at model with index 0 */
  pll_set_subst_params(partition, 0, 0, &(model.substitution_rates()[0]));

  /* set rate categories */
  pll_set_category_rates(partition, rate_cats);

  return make_tuple(partition, tree);

}
