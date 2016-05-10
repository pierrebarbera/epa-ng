#pragma once

#include <memory>

#include "pllhead.hpp"
#include "Sequence.hpp"
#include "constants.hpp"
#include "Model.hpp"
#include "Range.hpp"
#include "Placement.hpp"
#include "Tree.hpp"
#include "pll_util.hpp"

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
class Tiny_Tree
{
public:
  Tiny_Tree(pll_utree_t * edge_node, unsigned int branch_id, Tree& reference_tree, Model model,
    bool opt_branches) : Tiny_Tree(edge_node, branch_id, reference_tree,
    model, opt_branches, Range(0, reference_tree.partition()->sites), false) {};
  Tiny_Tree(pll_utree_t * edge_node , unsigned int branch_id, Tree& reference_tree, Model model,
    bool opt_branches, Range reference_tip_range, bool ranged);

  ~Tiny_Tree() = default;

  Tiny_Tree(Tiny_Tree const& other) = delete;
  Tiny_Tree(Tiny_Tree&& other) = default;

  Tiny_Tree& operator= (Tiny_Tree const& other) = delete;
  Tiny_Tree& operator= (Tiny_Tree && other) = default;

  Placement place(const Sequence& s);

  void opt_branches(bool b) {opt_branches_ = b;};

private:
  // pll structures
  std::unique_ptr<pll_partition_t, partition_deleter> partition_;
  std::unique_ptr<pll_utree_t, utree_deleter> tree_;

  bool opt_branches_;
  double original_branch_length_;
  Model model_;
  Range reference_tip_range_;
  bool tip_tip_case_ = false;
  bool ranged_computation_;
  unsigned int branch_id_;

};
