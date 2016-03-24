#pragma once

#include "pllhead.hpp"
#include "Sequence.hpp"
#include "constants.hpp"
#include "Model.hpp"
#include "Range.hpp"
#include "Placement.hpp"

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
  Tiny_Tree(pll_utree_t * edge_node, unsigned int branch_id, pll_partition_t * old_partition, Model model,
    bool opt_branches) : Tiny_Tree(edge_node, branch_id, old_partition,
    model, opt_branches, Range(0, old_partition->sites), false) {};
  Tiny_Tree(pll_utree_t * edge_node , unsigned int branch_id, pll_partition_t * old_partition, Model model,
    bool opt_branches, Range reference_tip_range, bool ranged);

  Tiny_Tree(Tiny_Tree const& other) = delete;

  Tiny_Tree(Tiny_Tree&& other)
    : partition_(other.partition_), tree_(other.tree_),
      opt_branches_(other.opt_branches_), original_branch_length_(other.original_branch_length_),
      model_(other.model_),
      reference_tip_range_(other.reference_tip_range_), tip_tip_case_(other.tip_tip_case_),
      ranged_computation_(other.ranged_computation_)
  {
    other.partition_ = nullptr;
    other.tree_ = nullptr;
  }
  ~Tiny_Tree();

  Tiny_Tree& operator = (Tiny_Tree const& other) = delete;
  Tiny_Tree& operator = (Tiny_Tree && other) = delete;

  // returns, in order, : likelihood, distal length, pendant length
  Placement place(const Sequence& s);

  void opt_branches(bool b) {opt_branches_ = b;};

private:
  // pll structures
  pll_partition_t * partition_ = nullptr;
  pll_utree_t * tree_ = nullptr;

  bool opt_branches_;
  const double original_branch_length_;
  Model model_;
  Range reference_tip_range_;
  bool tip_tip_case_ = false;
  bool ranged_computation_;
  unsigned int branch_id_;

};
