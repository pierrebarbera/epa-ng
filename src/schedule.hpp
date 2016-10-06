#pragma once

#include <vector>
#include <stdexcept>
#include <cassert>
#include <cmath>
#include <numeric>
#include <algorithm>

void to_difficulty(std::vector<double> perstage_avg)
{
  auto min = std::min_element(perstage_avg.begin(), perstage_avg.end());
  for_each(perstage_avg.begin(), perstage_avg.end(), [min](double& x){x /= min};);
}

std::vector<unsigned int> solve(unsigned int stages, unsigned int nodes, std::vector<double> difficulty_per_stage)
{
  assert(difficulty_per_stage.size() == stages);
  if (nodes < stages)
    throw std::runtime_error{"Must have more or equal number of nodes than stages"};

  std::vector<unsigned int> nodes_per_stage(stages);
  
  auto x1 = std::accumulate(difficulty_per_stage.begin(), difficulty_per_stage.end(), 0.0);
  x1 = nodes / x1;

  for (unsigned int i = 0; i < stages; ++i)
  {
    nodes_per_stage[i] = ceil(difficulty_per_stage[i] * x1);
  }

  int off_by;

  while ( (off_by = std::accumulate(nodes_per_stage.begin(), nodes_per_stage.end(), 0) - nodes) ) 
  {
    auto max_stage = std::max_element(nodes_per_stage.begin(), nodes_per_stage.end());
    if (off_by < 0)
    {
      *max_stage += 1;
    }
    else 
    {
      *max_stage -= 1;
    }
  }

  return nodes_per_stage;
}