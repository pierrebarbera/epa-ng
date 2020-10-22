#include "core/BranchBuffer.hpp"

#include "core/pll/epa_pll_util.hpp"
#include "core/Work.hpp"

#include <algorithm>
#include <future>

/** calculate the next block of branches, according to the branch traversal in
 * in the ref tree
 */
static void calc_block( BranchBuffer::container_type& buffer,
                        size_t const requested_size,
                        Tree* tree,
                        size_t* traversal_start,
                        BranchBuffer::whitelist_type const& treat_branch )
{
  bool const use_whitelist = not treat_branch.empty();

  auto& memsave = tree->memsave();
  assert( memsave );
  auto const& traversal = memsave.traversal();
  auto& branch_id       = tree->branch_id();

  // how many branches are there still in the complete traversal, until we've
  // calculated them all (in this overall iteration)
  size_t const branches_left = traversal.size() - *traversal_start;

  size_t const block_size = std::min( branches_left, requested_size );
  buffer.clear();
  buffer.reserve( block_size );

  for( size_t i = *traversal_start; i < *traversal_start + block_size; ++i ) {
    // for the current branch, according to the traversal
    auto branch_node = traversal[ i ];
    auto const cur_branch_id  = branch_id[ branch_node->node_index ];
    // skip the branch if it isn't whitelisted
    if( use_whitelist and not treat_branch[ cur_branch_id ] ) {
      continue;
    }
    // compute the required CLVs
    partial_compute_clvs( tree->tree(),
                          tree->nums(),
                          memsave.subtree_sizes(),
                          branch_node,
                          tree->partition() );

    // and persist them into a dedicated object holding a copy of the CLV
    // buffers for a given branch
    buffer.emplace_back( branch_node,
                         cur_branch_id,
                         *tree,
                         true ); // deep-copy CLVs
  }

  // update the start of the next traversal to be one beyond the block we just
  // finished
  *traversal_start += block_size;
}

/**
 * Whitelisted BranchBuffer: skip branches that don't occur in given Work object
 */
BranchBuffer::BranchBuffer( Tree* tree,
                            size_t const block_size,
                            Work const& work )
    : tree_( tree )
    , buffer_size_( block_size )
    , traversal_start_( 0 )
    , whitelist_( work.empty() ? 0ul : tree->nums().branches, false )
{
  if( not work.empty() ) {
    for( auto iter = work.bin_cbegin(); iter != work.bin_cend(); ++iter)
    {
      auto const branch_id    = iter->first;
      whitelist_[ branch_id ] = true;
    }
  }

  prefetcher_ = std::async( std::launch::async,
                            calc_block,
                            std::ref( buffer_ ),
                            buffer_size_,
                            tree_,
                            &traversal_start_,
                            std::ref( whitelist_ ) );
}

/**
 * Create the object and immediately start pre-calculating the first block
 */
BranchBuffer::BranchBuffer( Tree* tree, size_t const block_size )
    : BranchBuffer( tree, block_size, Work() )

{
}

size_t BranchBuffer::get_next( BranchBuffer::container_type& result )
{
  // join prefetching thread to ensure new block exists
  if( prefetcher_.valid() ) {
    prefetcher_.wait();
  }

  // perform pointer swap to data
  std::swap( result, buffer_ );

  // kick off calculation of the next block
  prefetcher_ = std::async( std::launch::async,
                            calc_block,
                            std::ref( buffer_ ),
                            buffer_size_,
                            tree_,
                            &traversal_start_,
                            std::ref( whitelist_ ) );

  // return the size of the returned block
  return result.size();
}

BranchBuffer::~BranchBuffer()
{
  // avoid dangling threads
  if( prefetcher_.valid() ) {
    prefetcher_.wait();
  }
}
