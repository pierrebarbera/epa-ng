#ifndef EPA_TINY_TREE_H_
#define EPA_TINY_TREE_H_

#include "pllhead.hpp"
#include "Sequence.hpp"

constexpr unsigned int OLD_LEFT_CLV_INDEX = 0;
constexpr unsigned int OLD_RIGHT_CLV_INDEX = 1;
constexpr unsigned int NEW_TIP_CLV_INDEX = 2;
constexpr unsigned int INNER_CLV_INDEX = 3;

constexpr unsigned int NUM_OPS = 1;

/* Encapsulates a smallest possible unrooted tree (3 tip nodes, 1 inner node)
  for use in edge insertion:

             [2]
           new_tip
              |
              |
            inner [3]
           /     \
          /       \
      old_left  old_right
        [0]       [1]

  where old left/right are the nodes adjacent to the insertion edge in the reference tree.
  new_tip represents the newly added sequence.

*/
class Tiny_Tree
{
public:
  Tiny_Tree(pll_utree_t * edge_node, pll_partition_t * old_partition, bool heuristic=true);

  Tiny_Tree(Tiny_Tree const& other) = delete;

  Tiny_Tree(Tiny_Tree&& other)
    : partition_(other.partition_)
    , ops_(other.ops_)
    , old_branch_length_ (other.old_branch_length_)
    , heuristic_(other.heuristic_)
  {
    other.partition_ = nullptr;
  }

  Tiny_Tree& operator = (Tiny_Tree const& other) = delete;
  Tiny_Tree& operator = (Tiny_Tree && other) = delete;

  ~Tiny_Tree();

  double place(const Sequence& s);

private:
  // pll structures
  pll_partition_t * partition_ = nullptr;
  // pll_utree_t * tree_;
  pll_operation_t ops_;

  // useful numbers
  double old_branch_length_;
  const bool heuristic_;

};

#endif
