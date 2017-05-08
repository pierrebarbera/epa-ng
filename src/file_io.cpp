#include "file_io.hpp"

#include <stdexcept>
#include <fstream>

#include "pllhead.hpp"
#include "MSA.hpp"
#include "constants.hpp"
#include "pll_util.hpp"

using namespace std;

/* reads in sequences from a file
  msa_file: string specifying the file path
*/
MSA build_MSA_from_file(const string& msa_file)
{
  /* open the file */
  auto file = pll_fasta_open(msa_file.c_str(), pll_map_fasta);
  if (!file) {
    throw runtime_error{string("Cannot open file ") + msa_file};
  }

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

  free(sequence);
  free(header);

  /* read the rest */
  while (pll_fasta_getnext(file, &header, &header_length, &sequence, &sequence_length, &sequence_number))
  {

    if (sites && sites != sequence_length)
      throw runtime_error{"MSA file does not contain equal size sequences"};

    if (!sites) sites = sequence_length;

    msa.append(header, sequence);
    free(sequence);
    free(header);
  }

  if (pll_errno != PLL_ERROR_FILE_EOF)
    throw runtime_error{string("Error while reading file: ") +  msa_file};

  pll_fasta_close(file);

  return msa;
}

pll_utree_t * build_tree_from_file(const string& tree_file, Tree_Numbers& nums)
{
  unsigned int num_tip_nodes;

  pll_utree_t * tree;
  pll_rtree_t * rtree;

  // load the tree unrooted
  if (!(rtree = pll_rtree_parse_newick(tree_file.c_str(), &num_tip_nodes))) {
   if (!(tree = pll_utree_parse_newick(tree_file.c_str(), &num_tip_nodes))) {
     throw runtime_error{"Treeparsing failed!"};
   }
  } else {
   tree = pll_rtree_unroot(rtree);
   pll_rtree_destroy(rtree, nullptr);

   /* optional step if using default PLL clv/pmatrix index assignments */
   pll_utree_reset_template_indices(tree, num_tip_nodes);
  }

  if (num_tip_nodes < 3) {
    throw runtime_error{"Number of tip nodes too small"};
  }

  nums = Tree_Numbers(num_tip_nodes);

  set_missing_branch_lengths(tree, DEFAULT_BRANCH_LENGTH);

  return tree;
}

pll_partition_t *  build_partition_from_file( const Model& model, 
                                              Tree_Numbers& nums, 
                                              const int num_sites)
{
  assert(nums.tip_nodes); // nums must have been initialized correctly

  auto attributes = PLL_ATTRIB_ARCH_CPU;
#ifdef __AVX
  attributes = PLL_ATTRIB_ARCH_AVX;
#elif __SSE3
  attributes = PLL_ATTRIB_ARCH_SSE;
#endif
  // attributes |= PLL_ATTRIB_PATTERN_TIP;
  attributes |= PLL_ATTRIB_SITES_REPEATS;

  auto partition = pll_partition_create(nums.tip_nodes,
           nums.inner_nodes * 3, //number of extra clv buffers: 3 for every direction on the node
           model.states(),
           num_sites,
           1,
           nums.branches,
           model.rate_cats(),
           (nums.inner_nodes * 3) + nums.tip_nodes, /* number of scaler buffers */
           attributes);

  if (!partition) {
    throw runtime_error{std::string("Could not create partition (build_partition_from_file). pll_errmsg: ") + pll_errmsg};
  }

  vector<double> rate_cats(model.rate_cats(), 0.0);

  /* compute the discretized category rates from a gamma distribution
     with alpha shape */
  pll_compute_gamma_cats(model.alpha(), model.rate_cats(), &rate_cats[0]);
  pll_set_frequencies(partition, 0, &(model.base_frequencies()[0]));
  pll_set_subst_params(partition, 0, &(model.substitution_rates()[0]));
  pll_set_category_rates(partition, &rate_cats[0]);

  return partition;

}

void file_check(const string& file_path)
{
  ifstream file(file_path.c_str());
  if (!file.good()) {
    throw runtime_error{string("file_check failed: ") + file_path};
  }

  file.close();
}
