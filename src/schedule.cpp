#include "schedule.hpp"

#include <stdexcept>
#include <cassert>
#include <cmath>
#include <numeric>
#include <algorithm>

void to_difficulty(std::vector<double>& perstage_avg)
{
  auto min = *std::min_element(perstage_avg.begin(), perstage_avg.end());
  for_each(perstage_avg.begin(), perstage_avg.end(), 
    [min](double& x){x /= min;}
    );
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

void assign(std::vector<unsigned int>& nodes_per_stage, 
            std::unordered_map<int, std::unordered_map<int, int>>& rank_assignm,
            int* local_stage)
{
  int rank = 0;
  for (unsigned int stage = 0; stage < nodes_per_stage.size(); ++stage)
  {
    auto nodes = nodes_per_stage[stage];
    for (unsigned int j = 0; j < nodes; ++j)
    {
      rank_assignm[stage][j] = rank;
      rank++;
    }
  }
}