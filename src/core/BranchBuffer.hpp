#pragma once

#include "core/Work.hpp"
#include "core/pll/pllhead.hpp"
#include "tree/Tiny_Tree.hpp"
#include "tree/Tree.hpp"

#include <future>
#include <unordered_map>
#include <vector>

/**
 * Buffers x branches / TinyTrees as a block for computation, while
 * pre-calculating the next block of x, asynchronously.
 */
class BranchBuffer {
  public:
  using container_type = std::vector< Tiny_Tree >;
  using whitelist_type = std::vector< bool >;

  BranchBuffer() = delete;

  BranchBuffer( Tree* tree, size_t const block_size );
  BranchBuffer( Tree* tree, size_t const block_size, Work const& work );

  // allow moving, but not copying
  BranchBuffer( BranchBuffer&& other ) = default;
  BranchBuffer& operator=( BranchBuffer&& other ) = default;

  BranchBuffer( BranchBuffer const& other ) = delete;
  BranchBuffer& operator=( BranchBuffer const& other ) = delete;

  // explicitly defined in .cpp, waits for dangling async
  ~BranchBuffer();

  size_t get_next( container_type& result );

  Tree const& tree() { return *tree_; }
  whitelist_type const& whitelist() const { return whitelist_; }

  private:
  Tree* tree_ = nullptr;
  size_t buffer_size_;
  container_type buffer_;
  #ifdef __PREFETCH
  std::future< void > prefetcher_;
  #endif
  size_t traversal_start_;
  whitelist_type whitelist_;
};
