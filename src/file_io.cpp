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
  if (!file)
    throw runtime_error{string("Cannot open file ") + msa_file};

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

  auto tree = pll_utree_parse_newick(tree_file.c_str(), &num_tip_nodes);

  if (num_tip_nodes < 3)
    throw runtime_error{"Number of tip nodes too small"};

  nums.init(num_tip_nodes);

  set_missing_branch_lengths(tree, DEFAULT_BRANCH_LENGTH);

  return tree;
}

pll_partition_t *  build_partition_from_file(const Model& model, Tree_Numbers& nums,  const int num_sites)
{
  assert(nums.tip_nodes); // nums must have been initialized correctly

  auto attributes = PLL_ATTRIB_ARCH_SSE;
#ifdef __AVX
  attributes = PLL_ATTRIB_ARCH_AVX;
#endif

  attributes |= PLL_ATTRIB_PATTERN_TIP;

  auto partition = pll_partition_create(nums.tip_nodes,
           nums.inner_nodes * 3, //number of extra clv buffers: 3 for every direction on the node
           STATES,
           num_sites,
           0, // number of mixture models
           1,
           nums.branches,
           RATE_CATS,
           (nums.inner_nodes * 3) + nums.tip_nodes, /* number of scaler buffers */
           pll_map_nt,
           attributes);

  double rate_cats[RATE_CATS] = {0};

  /* compute the discretized category rates from a gamma distribution
     with alpha shape */
  pll_compute_gamma_cats(model.alpha(), RATE_CATS, rate_cats);
  pll_set_frequencies(partition, 0, 0, &(model.base_frequencies()[0]));
  pll_set_subst_params(partition, 0, 0, &(model.substitution_rates()[0]));
  pll_set_category_rates(partition, rate_cats);

  return partition;

}

// static pll_partition_t * create_out_of_core_partition(FILE * bin_file)
// {
//   auto attributes = PLL_BINARY_ATTRIB_PARTITION_DUMP_CLV | PLL_BINARY_ATTRIB_PARTITION_DUMP_WGT;
//   auto partition = pll_binary_partition_load(bin_file, nullptr, attributes, pll_map_nt);
//   if (!partition)
//     throw runtime_error{"Error while building partition from binary file."};
//
//   return partition;
// }
//
// pll_partition_t *  build_partition_from_binary(const string& bin_file, const bool out_of_core)
// {
//   // savely open the file
//   FILE * bin_file_ptr;
//   if(!(bin_file_ptr = fopen(bin_file.c_str(), "r")))
//     throw runtime_error{"Couldn't open binary file."};
//
//   // decide whether to load the full partition, clv's and all
//   pll_partition_t * shell_partition = out_of_core ? create_out_of_core_partition(bin_file_ptr) : nullptr;
//
//   auto attributes = PLL_BINARY_ATTRIB_PARTITION_DUMP_CLV | PLL_BINARY_ATTRIB_PARTITION_DUMP_WGT;
//   if (out_of_core)
//     attributes = 0;
//
//   auto partition = pll_binary_partition_load(bin_file_ptr, shell_partition, attributes, pll_map_nt);
//
//   fclose(bin_file_ptr);
//
//   return partition;
// }

void file_check(const string& file_path)
{
  ifstream file(file_path.c_str());
  if (!file.good())
    throw runtime_error{string("file_check failed: ") + file_path};

  file.close();
}
