#include "io/file_io.hpp"

#include <stdexcept>
#include <fstream>
#include <functional>

#include "core/pll/pllhead.hpp"
#include "core/pll/pll_util.hpp"
#include "core/pll/rtree_mapper.hpp"
#include "core/raxml/Model.hpp"
#include "io/msa_reader.hpp"
#include "seq/MSA.hpp"
#include "util/constants.hpp"
#include "util/logging.hpp"


typedef struct fasta_record_s {
  pll_fasta_t* file;
  char * sequence = NULL;
  char * header = NULL;
} fasta_record_t;

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

MSA build_MSA_from_file(const std::string& msa_file,
                        const MSA_Info& info,
                        const bool premasking)
{
  MSA msa;
  auto reader = make_msa_reader(msa_file, info, premasking, false);
  reader->read_next(msa, std::numeric_limits<size_t>::max());

  return msa;
}

static void recurse_post_order(pll_unode_t const * const node,
                        std::vector<unsigned int>& translation,
                        unsigned int * const rooted_index)
{
  if ( not node->next ) {
    translation.push_back( *rooted_index );
    *rooted_index = *rooted_index + 1;
    return;
  }

  recurse_post_order( node->next->back,       translation, rooted_index );
  recurse_post_order( node->next->next->back, translation, rooted_index );

  translation.push_back( *rooted_index );
  *rooted_index = *rooted_index + 1;
}

static rtree_mapper determine_edge_num_translation(pll_unode_t const * const vroot,
                                            bool const left,
                                            double const proximal_length,
                                            double const distal_length)
{
  // map from unrooted edge nums to rooted edge nums
  rtree_mapper::map_type unrooted_to_rooted;
  unsigned int rooted_index = 0;

  unsigned int rtree_proximal_edge  = 0;
  unsigned int rtree_distal_edge = 0;
  unsigned int utree_root_edge  = 0;

  /* since we are in a special situation of a) starting at a top level trifurcation
   * and b) having to treat one edge in a special way, we do part of the normally
   * recursive legwork here already
   */

  if ( left ) {
    // if the vroot sits on the left subtree of the original rtree, we can recurse twice,
    // then account for the missing branch w.r.t. the rooted tree
    recurse_post_order( vroot->back,        unrooted_to_rooted, &rooted_index );
    recurse_post_order( vroot->next->back,  unrooted_to_rooted, &rooted_index );

    // remember the edge num in the rtree of the branch leading to left subtree
    rtree_proximal_edge = rooted_index;
    rooted_index++;

    // keep recursing down the originally right path
    recurse_post_order( vroot->next->next->back, unrooted_to_rooted, &rooted_index );

    // remember the edge num in the rtree of the branch leading to the right subtree
    rtree_distal_edge = unrooted_to_rooted.back();

    // in this case the edge num in the utree that previously held the rtree root is the last branch by definiton
    utree_root_edge = unrooted_to_rooted.size() - 1u;
  } else {
    // if the vroot sits on the right subtree, we can recurse three times normally, then in the
    // end track the extraneous edge num

    // remember the edge num in the rtree of the branch leading to left subtree
    rtree_distal_edge = rooted_index;
    // ... which is also the utree edge num of the branch that held the rtree root
    utree_root_edge = rtree_distal_edge;

    recurse_post_order( vroot->back,  unrooted_to_rooted, &rooted_index );
    recurse_post_order( vroot, unrooted_to_rooted, &rooted_index );

    // remember the edge num in the rtree of the branch leading to right subtree
    // this also removes the uneeded mapping for the edge that was traversed twice
    // (the edge where the rtree root was)
    rtree_proximal_edge = unrooted_to_rooted.back();
    unrooted_to_rooted.pop_back();

  }

  rtree_mapper mapper(utree_root_edge,
                      rtree_proximal_edge,
                      rtree_distal_edge,
                      proximal_length,
                      distal_length,
                      left );
  mapper.map( std::move(unrooted_to_rooted) );
  return mapper;
}

pll_utree_s * build_tree_from_file( const std::string& tree_file,
                                    Tree_Numbers& nums,
                                    rtree_mapper& mapper)
{
  pll_utree_t * tree;
  pll_rtree_t * rtree;

  // load the tree unrooted
  // if we can load it unrooted, then do so
  if ( (tree = pll_utree_parse_newick(tree_file.c_str()))) {

  // otherwise try to parse rooted, and unroot the tree
  } else if ( ( rtree = pll_rtree_parse_newick( tree_file.c_str() ) ) ) {

    // convert the tree
    tree = pll_rtree_unroot( rtree );

    // is the virtual root on the left side of the rtree?
    // (this is the case if the left child node isn't a tip)
    bool const left = rtree->root->left->left;

    // what was the original branch length that used to be attached to the now vroot?
    // lets by default assume the "utree distal" edge is the right one
    double distal_length  = rtree->root->left->length;
    double proximal_length = rtree->root->right->length;

    if ( left ) {
      /* if the virtual root sits on the left subtree we need to correct for an inconistent way the pll_rtee_unroot
       * function returns the root node. It returns the node whose ->back is the right hand subtree (from the
       * rooted perspective), which would be wrong as this is the subtree that post and pre order
       * traversals recurse into first.
       */
      tree->vroot = tree->vroot->next;

      // update branch lengths: they are flipped in this case
      double tmp = distal_length;
      distal_length = proximal_length;
      proximal_length = tmp;
    }

    // get the translation
    mapper = determine_edge_num_translation(  tree->vroot,
                                              left,
                                              proximal_length,
                                              distal_length );

    pll_rtree_destroy(rtree, nullptr);


    /* optional step if using default PLL clv/pmatrix index assignments */
    pll_utree_reset_template_indices( tree->vroot, tree->tip_count );

  // if both fails, abort
  } else {
    throw std::runtime_error{std::string("Treeparsing failed! ") + pll_errmsg};
  }

  if (not tree->binary) {
    throw std::invalid_argument{"Input Tree contains multifurcations (polytomies)!"};
  }

  if (tree->tip_count < 3) {
    throw std::runtime_error{"Number of tip nodes too small"};
  }

  nums = Tree_Numbers(tree->tip_count);

  set_missing_branch_lengths(tree, DEFAULT_BRANCH_LENGTH);

  return tree;
}

static unsigned int simd_autodetect()
{
  if (PLL_STAT(avx2_present))
    return PLL_ATTRIB_ARCH_AVX2;
  else if (PLL_STAT(avx_present))
    return PLL_ATTRIB_ARCH_AVX;
  else if (PLL_STAT(sse3_present))
    return PLL_ATTRIB_ARCH_SSE;
  else
    return PLL_ATTRIB_ARCH_CPU;
}

pll_partition_t *  make_partition(const raxml::Model& model,
                                  Tree_Numbers& nums,
                                  const int num_sites,
                                  const Options options)
{
  assert(nums.tip_nodes); // nums must have been initialized correctly

  auto attributes = simd_autodetect();

  if ( (options.scaling == Options::NumericalScaling::kOn) or
     ( (options.scaling == Options::NumericalScaling::kAuto) and nums.large_tree() ) ) {
    attributes = PLL_ATTRIB_RATE_SCALERS;
  }

  if ( options.repeats ) {
    attributes |= PLL_ATTRIB_SITE_REPEATS;
  } else {
    attributes |= PLL_ATTRIB_PATTERN_TIP;
  }

  auto partition = pll_partition_create(nums.tip_nodes,
           nums.inner_nodes * 3, //number of extra clv buffers: 3 for every direction on the node
           model.num_states(),
           num_sites,
           1,
           nums.branches,
           model.num_ratecats(),
           (nums.inner_nodes * 3) + nums.tip_nodes, /* number of scaler buffers */
           attributes);

  if (not partition) {
    throw std::runtime_error{std::string("Could not create partition (make_partition). pll_errmsg: ") + pll_errmsg};
  }

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
