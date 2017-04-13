#include "Binary.hpp"

#include <stdexcept>
#include <algorithm>

#include "constants.hpp"
#include "Tree.hpp"

using namespace std;

int safe_fclose(FILE* fptr) { return fptr ? fclose(fptr) : 0; }

Binary::Binary(Binary && other) : bin_fptr_(nullptr, safe_fclose)
{
  bin_fptr_ = std::move(other.bin_fptr_);
  map_ = std::move(other.map_);
}

Binary& Binary::operator=(Binary && other)
{
  bin_fptr_ = std::move(other.bin_fptr_);
  map_ = std::move(other.map_);
  return *this;
}

Binary::Binary(const string& binary_file_path) : bin_fptr_(nullptr, safe_fclose)
{
  // open the binary file
  pll_binary_header_t header;
  bin_fptr_ = file_ptr_type(pllmod_binary_open(binary_file_path.c_str(), &header), safe_fclose);

  if (!bin_fptr_)
    throw runtime_error{"Could not open binary file for reading."};

  if (header.access_type != PLLMOD_BIN_ACCESS_RANDOM)
    throw runtime_error{"Binary file must be random access enabled."};

  if (header.n_blocks <= 0)
    throw runtime_error{string("Binary file header must have nonzero positive number of blocks: ")
              + to_string(header.n_blocks)};

  // proccess the random access map
  unsigned int n_blocks;
  pll_block_map_t* block_map;
  block_map = pllmod_binary_get_map(bin_fptr_.get(), &n_blocks);

  assert(block_map);
  assert(n_blocks);

  for (size_t i = 0; i < n_blocks; i++)
    map_.push_back(block_map[i]);

  free(block_map);
}

static long int get_offset(vector<pll_block_map_t>& map, const int block_id)
{
  auto item = map.begin();
  while (item != map.end())
  {
    if (item->block_id == block_id) {
      break;
    }
    item++;
  }

  if(item == map.end())
    throw runtime_error{string("Map does not contain block_id: ") + to_string(block_id)};
  return item->block_offset;
}

void Binary::load_clv(pll_partition_t * partition, const unsigned int clv_index)
{
  std::lock_guard<std::mutex> lock(file_mutex_);
  assert(bin_fptr_);
  assert(clv_index < partition->clv_buffers + partition->tips);
  if (partition->attributes & PLL_ATTRIB_PATTERN_TIP)
    assert(clv_index >= partition->tips);

  if (!(partition->clv[clv_index]))
  {
    size_t clv_size = partition->sites * partition->states_padded *
                      partition->rate_cats*sizeof(double);
    partition->clv[clv_index] = (double*) pll_aligned_alloc(clv_size, partition->alignment);
    if (!partition->clv[clv_index])
      throw runtime_error{"Could not allocate CLV memory"};
  }

  unsigned int attributes;
  auto err = pllmod_binary_clv_load(
    bin_fptr_.get(),
    0,
    partition,
    clv_index,
    &attributes,
    get_offset(map_, clv_index));

  if (err != PLL_SUCCESS)
    throw runtime_error{string("Loading CLV failed: ") + pll_errmsg + std::string(". CLV index: ") + std::to_string(clv_index)};

}

void Binary::load_tipchars(pll_partition_t * partition, const unsigned int tipchars_index)
{
  std::lock_guard<std::mutex> lock(file_mutex_);
  assert(bin_fptr_);
  assert(tipchars_index < partition->tips);
  assert(partition->attributes & PLL_ATTRIB_PATTERN_TIP);

  unsigned int type, attributes;
  size_t size;

  auto ptr = pllmod_binary_custom_load(bin_fptr_.get(), 0, &size, &type, &attributes, get_offset(map_, tipchars_index));
  if (!ptr)
    throw runtime_error{string("Loading tipchar failed: ") + pll_errmsg};

  partition->tipchars[tipchars_index] = (unsigned char*)ptr;
}

void Binary::load_scaler(pll_partition_t * partition, const unsigned int scaler_index)
{
  std::lock_guard<std::mutex> lock(file_mutex_);
  assert(bin_fptr_);
  assert(scaler_index < partition->scale_buffers);

  auto block_offset = partition->clv_buffers + partition->tips;

  unsigned int type, attributes;
  size_t size;

  auto ptr = pllmod_binary_custom_load(bin_fptr_.get(), 0,
    &size, &type, &attributes, get_offset(map_, block_offset + scaler_index));
  if (!ptr)
    throw runtime_error{string("Loading scaler failed: ") + pll_errmsg};

  partition->scale_buffer[scaler_index] = (unsigned int*)ptr;
}

pll_partition_t* Binary::load_partition()
{
  std::lock_guard<std::mutex> lock(file_mutex_);
  // make skeleton partition that only allocates the pointers to the clv/tipchar buffers
  unsigned int attributes = PLLMOD_BIN_ATTRIB_PARTITION_LOAD_SKELETON;
  auto partition =  pllmod_binary_partition_load(bin_fptr_.get(), 0, nullptr, &attributes,
    get_offset(map_, -1));

  if (!partition)
    throw runtime_error{string("Error loading partition: ") + pll_errmsg};

  return partition;
}

pll_utree_t* Binary::load_utree()
{
  std::lock_guard<std::mutex> lock(file_mutex_);
  unsigned int attributes = 0;
  auto tree =  pllmod_binary_utree_load(bin_fptr_.get(), 0, &attributes, get_offset(map_, -2));
  if (!tree)
    throw runtime_error{string("Loading tree: ") + pll_errmsg};
  return tree;
}

/**
  Writes the structures and data encapsulated in Tree to the specified file in the binary format.
  Writes them in such a way that the Binary class can read them.
*/
void dump_to_binary(Tree& tree, const string& file)
{
  auto num_clvs = tree.partition()->clv_buffers;
  auto num_tips = tree.partition()->tips;
  auto num_scalers = tree.partition()->scale_buffers;
  auto max_clv_index = num_clvs + num_tips;

  pll_binary_header_t header;
  auto fptr =  pllmod_binary_create(
    file.c_str(),
    &header,
    PLLMOD_BIN_ACCESS_RANDOM,
    2 + num_clvs + num_tips + num_scalers);

  if(!fptr)
    throw runtime_error{string("Opening binary file for writing: ") + pll_errmsg};

  unsigned int attributes = PLLMOD_BIN_ATTRIB_UPDATE_MAP | PLLMOD_BIN_ATTRIB_PARTITION_DUMP_WGT;

  int block_id = -2;
  bool use_tipchars = tree.partition()->attributes & PLL_ATTRIB_PATTERN_TIP;

  // dump the utree structure
  if(!pllmod_binary_utree_dump(fptr, block_id++, tree.tree(), num_tips, attributes))
    throw runtime_error{string("Dumping the utree to binary: ") + pll_errmsg};

  // dump the partition
  if(!pllmod_binary_partition_dump(fptr, block_id++, tree.partition(), attributes))
    throw runtime_error{string("Dumping partition to binary: ") + pll_errmsg};

  // dump the tipchars, but only if partition uses them
  size_t tip_index = 0;
  if (use_tipchars)
  {
    for (tip_index = 0; tip_index < num_tips; tip_index++)
    {
      if(!pllmod_binary_custom_dump(fptr, block_id++, tree.partition()->tipchars[tip_index],
      tree.partition()->sites * sizeof(unsigned char), attributes))
        throw runtime_error{string("Dumping tipchars to binary: ") + pll_errmsg};
    }
  }

  // dump the clvs
  for (size_t clv_index = tip_index; clv_index < max_clv_index; clv_index++)
  {
    if(!pllmod_binary_clv_dump(fptr, block_id++, tree.partition(), clv_index, attributes))
      throw runtime_error{string("Dumping clvs to binary: ") + pll_errmsg};
  }

  // dump the scalers
  for (size_t scaler_index = 0; scaler_index < num_scalers; scaler_index++)
  {
    if(!pllmod_binary_custom_dump(fptr, block_id++, tree.partition()->scale_buffer[scaler_index],
      tree.partition()->sites * sizeof(unsigned int), attributes))
      throw runtime_error{string("Dumping scalers to binary: ") + pll_errmsg};
  }

  fclose(fptr);
}
