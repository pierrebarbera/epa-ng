#pragma once

#include "core/pll/pllhead.hpp"
#include "tree/Tiny_Tree.hpp"
#include "tree/Tree.hpp"

#include <future>
#include <vector>

/**
 * Buffers x branches / TinyTrees as a block for computation, while pre-calculating
 * the next block of x, asynchronously.
 */
class BranchBuffer {
  public:
  using container_type = std::vector< TinyTree >;

  /**
   * Create the object and immediately start pre-calculating the first block
   */
  BranchBuffer( Tree* ref_tree,
                size_t const block_size );

  // allow moving, but not copying
  BranchBuffer( BranchBuffer&& other )      = default;
  BranchBuffer& operator=( BranchBuffer&& other )       = default;
  
  BranchBuffer( BranchBuffer const& other ) = delete;
  BranchBuffer& operator=( BranchBuffer const& other )  = delete;

  ~BranchBuffer();

  size_t get_next( container_type& result, size_t const block_size );

  private:
  Tree* ref_tree_ = nullptr;
  size_t buffer_size_;
  container_type buffer_;
  std::future< void > prefetcher_;
  size_t traversal_start_;
};
