#ifndef EPA_TREE_NUMBERS_H_
#define EPA_TREE_NUMBERS_H_

class Tree_Numbers
{
public:
  Tree_Numbers(){};
  ~Tree_Numbers() {};
  void init(int tn)
  {
    tip_nodes = tn;
    inner_nodes = tip_nodes - 2;
    nodes = inner_nodes + tip_nodes;
    branches = nodes - 1;
  };
  int tip_nodes;
  int inner_nodes;
  int nodes;
  int branches;
};

#endif
