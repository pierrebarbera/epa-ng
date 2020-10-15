#include "core/BranchBuffer.hpp"

#include "core/pll/epa_pll_util.hpp"

#include <algorithm>
#include <future>

/** calculate the next block of branches, according to the branch traversal in
 * in the ref tree
 */
static void calc_block( BranchBuffer::container_type& buffer,
                        size_t const requested_size,
                        Tree* tree,
                        size_t* traversal_start )
{
  auto& memsave = tree->memsave();
  assert( memsave );
  auto const& traversal = memsave.traversal();
  auto& branch_id       = tree->branch_id();

  // how many branches are there still in the complete traversal, until we've
  // calculated them all (in this overall iteration)
  size_t const branches_left = traversal.size() - *traversal_start;

  size_t const block_size = std::min( branches_left, requested_size );
  // buffer.clear();
  buffer.resize( block_size );

  for( size_t i = *traversal_start; i < *traversal_start + block_size; ++i ) {
    // for the current branch, according to the traversal
    auto branch_node = traversal[ i ];
    // compute the required CLVs
    partial_compute_clvs( tree->tree(),
                          tree->nums(),
                          memsave.subtree_sizes(),
                          branch_node,
                          tree->partition() );

    // and persist them into a dedicated object holding a copy of the CLV
    // buffers for a given branch
    size_t const buffer_index = i - *traversal_start;
    buffer[ buffer_index ]    = Tiny_Tree( branch_node,
                                        branch_id[ branch_node->node_index ],
                                        *tree,
                                        true ); // deep-copy CLVs
  }

  // update the start of the next traversal to be one beyond the block we just
  // finished
  *traversal_start += block_size;
}

/**
 * Create the object and immediately start pre-calculating the first block
 */
BranchBuffer::BranchBuffer( Tree* tree, size_t const block_size )
    : tree_( tree )
    , buffer_size_( block_size )
    , traversal_start_( 0 )
{
  prefetcher_ = std::async( std::launch::async,
                            calc_block,
                            std::ref( buffer_ ),
                            block_size,
                            tree_,
                            &traversal_start_ );
}

size_t BranchBuffer::get_next( BranchBuffer::container_type& result,
                               size_t const block_size )
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
                            block_size,
                            tree_,
                            &traversal_start_ );

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
