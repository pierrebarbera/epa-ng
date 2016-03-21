#pragma once

class Tree_Numbers
{
public:
  Tree_Numbers(){};
  ~Tree_Numbers() {};
  void init(unsigned int tn)
  {
    tip_nodes = tn;
    inner_nodes = tip_nodes - 2;
    nodes = inner_nodes + tip_nodes;
    branches = nodes - 1;
  };
  unsigned int tip_nodes = 0;
  unsigned int inner_nodes = 0;
  unsigned int nodes = 0;
  unsigned int branches = 0;
};
