#ifndef EPA_TINY_TREE_H_
#define EPA_TINY_TREE_H_

#include <tuple>

#include "pllhead.hpp"
#include "Sequence.hpp"
#include "constants.hpp"
#include "Model.hpp"

/* Encapsulates a smallest possible unrooted tree (3 tip nodes, 1 inner node)
  for use in edge insertion:

             [2]
           new_tip
              |
              |
            inner [3]
           /     \
          /       \
      proximal  distal
        [0]       [1]

  where proximal/distal are the nodes adjacent to the insertion edge in the reference tree.
  new_tip represents the newly added sequence.

*/
class Tiny_Tree
{
public:
  Tiny_Tree(pll_utree_t * edge_node, pll_partition_t * old_partition, Model model,
    bool opt_branches=false);

  Tiny_Tree(Tiny_Tree const& other) = delete;

  Tiny_Tree(Tiny_Tree&& other)
    : partition_(other.partition_), tree_(other.tree_), ops_(other.ops_),
      opt_branches_(other.opt_branches_)
  {
    other.partition_ = nullptr;
    other.tree_ = nullptr;
  }
  ~Tiny_Tree();

  Tiny_Tree& operator = (Tiny_Tree const& other) = delete;
  Tiny_Tree& operator = (Tiny_Tree && other) = delete;

  // returns, in order, : likelihood, distal length, pendant length
  std::tuple<double, double, double> place(const Sequence& s);

private:
  // pll structures
  pll_partition_t * partition_ = nullptr;
  pll_utree_t * tree_ = nullptr;
  pll_operation_t ops_;

  const bool opt_branches_;
  Model model_;

};

#endif
