#pragma once

#include <memory>
#include <utility>

#include "core/pll/pll_util.hpp"
#include "core/pll/pllhead.hpp"
#include "sample/Placement.hpp"
#include "seq/Sequence.hpp"
#include "tree/Tree.hpp"
#include "util/Options.hpp"
#include "util/constants.hpp"

/* Encapsulates a smallest possible unrooted tree (3 tip nodes, 1 inner node)
  for use in edge insertion:

             S:[-] C:[1]
           new_tip
              |
              |
            inner S:[1] C:[3]
           /     \
          /       \
      proximal    distal
    S:[0] C:[4]   S:[2] C:[2 or 5]

  where proximal/distal are the nodes adjacent to the insertion edge in the reference tree.
  new_tip represents the newly added sequence.

*/
class Tiny_Tree {
  public:
  Tiny_Tree( pll_unode_t* edge_node,
             unsigned int const branch_id,
             Tree& reference_tree,
             bool const opt_branches,
             Options const& options);

  Tiny_Tree()  = delete;
  ~Tiny_Tree() = default;

  Tiny_Tree( Tiny_Tree const& other ) = delete;
  Tiny_Tree( Tiny_Tree&& other )      = default;

  Tiny_Tree& operator=( Tiny_Tree const& other ) = delete;
  Tiny_Tree& operator=( Tiny_Tree&& other ) = default;

  Placement place( Sequence const& s );

  void get_persite_logl( char const nt, std::vector< double >& result ) const;

  // getters / setters

  unsigned int branch_id() const { return branch_id_; }
  double distal_length() const { return tree_->nodes[ 1 ]->length; }
  double proximal_length() const { return tree_->nodes[ 0 ]->length; }
  double pendant_length() const { return tree_->nodes[ 3 ]->length; }

  private:
  // pll structures
  std::unique_ptr< pll_partition_t, partition_deleter > partition_;
  std::unique_ptr< pll_utree_t, utree_deleter > tree_;

  bool opt_branches_;
  double original_branch_length_;
  bool premasking_ = true;
  bool sliding_blo_;
  unsigned int branch_id_;
};
