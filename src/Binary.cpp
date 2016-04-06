#include "Binary.hpp"

#include <stdexcept>

#include "constants.hpp"

using namespace std;

Binary::Binary(const string& binary_file_path) : bin_fptr_(nullptr)
{
  // open the binary file
  pll_binary_header_t header;
  bin_fptr_ = pll_binary_open(binary_file_path.c_str(), &header);

  if (!bin_fptr_)
    throw runtime_error{"Could not open binary file for reading."};

  if (header.access_type != PLL_BINARY_ACCESS_RANDOM)
    throw runtime_error{"Binary file must be random access enabled!"};

  // proccess the random access map
  // unsigned int n_blocks;
  // pll_block_map_t* block_map;
  // block_map = pll_binary_get_map(bin_fptr_, &n_blocks);
  //
  // for (size_t i = 0; i < n_blocks; i++)
  // {
  //   switch (block_map[i].block_id) {
  //     case PLL_BINARY_BLOCK_CLV:
  //       clv_offset_ = block_map[i].block_offset;
  //       break;
  //     case PLL_BINARY_BLOCK_PARTITION:
  //       partition_offset_ = block_map[i].block_offset;
  //       break;
  //     case PLL_BINARY_BLOCK_CUSTOM:
  //       tipchars_offset_ = block_map[i].block_offset;
  //       break;
  //     default:
  //       break;
  //   }
  // }
  //
  // // cleanup
  // free(block_map);
}

void Binary::load_clv(pll_partition_t * partition, const unsigned int clv_index)
{
  assert(bin_fptr_);
  assert(clv_index < partition->clv_buffers + partition->tips);
  assert(clv_index >= partition->tips);

  unsigned int attributes;
  auto err = pll_binary_clv_load(
    bin_fptr_,
    clv_index,
    partition,
    clv_index,
    &attributes,
    PLL_BINARY_ACCESS_SEEK);

  if (err != PLL_SUCCESS)
    throw runtime_error{string("Loading CLV failed with pll_errno: ") + to_string(pll_errno)};
}

void Binary::load_tipchars(pll_partition_t * partition, const unsigned int tipchars_index)
{
  assert(bin_fptr_);
  assert(tipchars_index < partition->tips);

  unsigned int type, attributes;
  size_t size;

  auto ptr = pll_binary_custom_load(bin_fptr_, tipchars_index, &size, &type, &attributes, PLL_BINARY_ACCESS_SEEK);
  if (!ptr)
    throw runtime_error{string("Loading tipchar failed with pll_errno: ") + to_string(pll_errno)};

  partition->tipchars[tipchars_index] = (char*)ptr;
}

void Binary::load_scaler(pll_partition_t * partition, const unsigned int scaler_index)
{
  assert(bin_fptr_);
  assert(scaler_index < partition->scale_buffers);

  auto block_offset = partition->clv_buffers + partition->tips;

  unsigned int type, attributes;
  size_t size;

  auto ptr = pll_binary_custom_load(bin_fptr_, block_offset + scaler_index, &size, &type, &attributes, PLL_BINARY_ACCESS_SEEK);
  if (!ptr)
    throw runtime_error{string("Loading scaler failed with pll_errno: ") + to_string(pll_errno)};

  partition->scale_buffer[scaler_index] = (unsigned int*)ptr;
}

static pll_partition_t* skeleton_partition()
{
  auto attributes = PLL_ATTRIB_ARCH_SSE;
#ifdef __AVX
  attributes = PLL_ATTRIB_ARCH_AVX;
#endif

  attributes |= PLL_ATTRIB_PATTERN_TIP;

  auto partition = pll_partition_create(
    0, // number of tip nodes
    0, // number of extra clv buffers
    STATES,
    0, // number of sites
    0, // number of mixture models
    1, // number of cincurrent subs. models
    0, // number of probabillity matrices
    RATE_CATS,
    0, // number of scale buffers
    pll_map_nt,
    attributes);

  if (!partition)
    throw runtime_error{string("Creating skeleton partition failed with pll_errno: ") + to_string(pll_errno)};

  // ensure clv, tipchar and scaler fields are only shallowly allocated


  return partition;
}

pll_partition_t* Binary::load_partition()
{
  // make skeleton partition that only allocates the pointers to the clv/tipchar buffers
  auto skelly = skeleton_partition();
  unsigned int attributes = PLL_BINARY_ATTRIB_PARTITION_DUMP_WGT;
  auto partition =  pll_binary_partition_load(bin_fptr_, -1, skelly, &attributes, pll_map_nt, 0);
  if (!partition)
    throw runtime_error{string("Loading partition failed with pll_errno: ") + to_string(pll_errno)};
  return partition;
}

pll_utree_t* Binary::load_utree()
{
  unsigned int attributes = 0;
  auto tree =  pll_binary_utree_load(bin_fptr_, -2, &attributes, 0);
  if (!tree)
    throw runtime_error{string("Loading tree failed with pll_errno: ") + to_string(pll_errno)};
  return tree;
}
