#include "file_io.hpp"

#include <stdexcept>
#include <fstream>

#include "pllhead.hpp"
#include "MSA.hpp"
#include "constants.hpp"
#include "pll_util.hpp"
#include "logging.hpp"


typedef struct fasta_record_s {
  pll_fasta_t* file;
  char * sequence = NULL;
  char * header = NULL;
} fasta_record_t;

static unsigned long length_till_newl(char* line)
{
  char c;
  unsigned long length = 0;
  while((c = *line++) and (c != '\n' or c != '\r')) {
    ++length;
  }
  return length;
}


int pll_fasta_fseek(pll_fasta_t* fd, const long int offset, const int whence)
{
  auto status = fseek(fd->fp, offset, whence);

  /* reset stripped char frequencies */
  fd->stripped_count = 0;
  for(size_t i=0; i<256; i++) {
    fd->stripped[i] = 0;
  }

  fd->line[0] = 0;
  if (!fgets(fd->line, PLL_LINEALLOC, fd->fp)) {
    pll_errno = PLL_ERROR_FILE_SEEK;
    snprintf(pll_errmsg, 200, "Unable to rewind and cache data");
    return -1;
  }
  fd->lineno = 1;

  return status;
}

static void epa_read_fasta( const std::string& msa_file, 
                            std::function<void(fasta_record_t&)> fn, 
                            const size_t offset=0, 
                            const size_t span=std::numeric_limits<size_t>::max())
{
  fasta_record_t record;
  /* open the file */
  record.file = pll_fasta_open(msa_file.c_str(), pll_map_fasta);
  if (!record.file) {
    throw std::runtime_error{std::string("Cannot open file ") + msa_file};
  }

  if(pll_fasta_fseek(record.file, offset, SEEK_SET)) {
    throw std::runtime_error{"Unable to fseek on the fasta file."};
  }

  LOG_DBG << record.file->line;

  long sequence_length;
  long header_length;
  long sequence_number;

  /* read sequences and make sure they are all of the same length */
  int sites = 0;
  size_t read_seqs = 0;

  /* read the first sequence seperately, so that the MSA object can be constructed */
  if (PLL_FAILURE == pll_fasta_getnext( record.file,
                                        &record.header,
                                        &header_length,
                                        &record.sequence,
                                        &sequence_length,
                                        &sequence_number)) {
    throw std::runtime_error{std::string("fasta_getnext failed: ") +  pll_errmsg};
  }

  sites = sequence_length;

  if (sites == -1 || sites == 0) {
    throw std::runtime_error{"Unable to read MSA record: invalid sequence length"};
  }

  fn(record);
  ++read_seqs;

  free(record.sequence);
  free(record.header);

  /* read the rest */
  while(read_seqs < span 
        and pll_fasta_getnext(record.file, 
                              &record.header, 
                              &header_length, 
                              &record.sequence, 
                              &sequence_length, 
                              &sequence_number))
  {

    if (sites && sites != sequence_length) {
      throw std::runtime_error{"MSA file does not contain equal size sequences"};
    }

    if (!sites) sites = sequence_length;

    fn(record);
    ++read_seqs;

    free(record.sequence);
    free(record.header);
  }

  // if (pll_errno != PLL_ERROR_FILE_EOF) {
  //   throw std::runtime_error{std::string("Error while reading file: ") +  msa_file};
  // }

  pll_fasta_close(record.file);

}

std::vector<size_t> get_offsets(const std::string& msa_file)
{
  std::vector<size_t> result;

  epa_read_fasta(msa_file, 
    [&](auto& record){
      result.push_back(ftell(record.file->fp)
        - length_till_newl(record.file->line));
  });

  return result;
}

MSA build_MSA_from_file(const std::string& msa_file, 
                        const size_t offset,
                        const size_t span)
{
  MSA msa;
  
  epa_read_fasta(
    msa_file, 
    [&](auto& record){
      msa.append(record.header, record.sequence);
    },
    offset, 
    span);

  return msa;
}

pll_utree_s * build_tree_from_file(const std::string& tree_file, Tree_Numbers& nums)
{
  pll_utree_t * tree;
  pll_rtree_t * rtree;

  // load the tree unrooted
  if (!(rtree = pll_rtree_parse_newick(tree_file.c_str()))) {
    if (!(tree = pll_utree_parse_newick(tree_file.c_str()))) {
      throw std::runtime_error{"Treeparsing failed!"};
    }
  } else {
    tree = pll_rtree_unroot(rtree);
    pll_rtree_destroy(rtree, nullptr);

    /* optional step if using default PLL clv/pmatrix index assignments */
    pll_utree_reset_template_indices(get_root(tree), tree->tip_count);
  }

  if (tree->tip_count < 3) {
    throw std::runtime_error{"Number of tip nodes too small"};
  }

  nums = Tree_Numbers(tree->tip_count);

  set_missing_branch_lengths(tree, DEFAULT_BRANCH_LENGTH);

  return tree;
}

pll_partition_t *  build_partition_from_file( const Model& model, 
                                              Tree_Numbers& nums, 
                                              const int num_sites,
                                              const bool repeats)
{
  assert(nums.tip_nodes); // nums must have been initialized correctly

  auto attributes = PLL_ATTRIB_ARCH_CPU;
#ifdef __AVX
  attributes = PLL_ATTRIB_ARCH_AVX;
#elif __SSE3
  attributes = PLL_ATTRIB_ARCH_SSE;
#endif

  if (repeats) {
    attributes |= PLL_ATTRIB_SITES_REPEATS;
  } else {
    attributes |= PLL_ATTRIB_PATTERN_TIP;
  }

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
    throw std::runtime_error{std::string("Could not create partition (build_partition_from_file). pll_errmsg: ") + pll_errmsg};
  }

  std::vector<double> rate_cats(model.rate_cats(), 0.0);

  /* compute the discretized category rates from a gamma distribution
     with alpha shape */
  pll_compute_gamma_cats( model.alpha(), 
                          model.rate_cats(), 
                          &rate_cats[0],
                          PLL_GAMMA_RATES_MEAN);
  pll_set_frequencies(partition, 
                      0, 
                      &(model.base_frequencies()[0]));
  pll_set_subst_params( partition, 
                        0, 
                        &(model.substitution_rates()[0]));
  pll_set_category_rates( partition, 
                          &rate_cats[0]);

  // if (repeats) {
  //   pll_resize_repeats_lookup(partition, ( REPEATS_LOOKUP_SIZE ) * 10);
  // }

  return partition;

}

void file_check(const std::string& file_path)
{
  std::ifstream file(file_path.c_str());
  if (!file.good()) {
    throw std::runtime_error{std::string("file_check failed: ") + file_path};
  }

  file.close();
}
