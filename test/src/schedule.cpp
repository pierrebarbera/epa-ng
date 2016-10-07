#include "Epatest.hpp"

#include "src/schedule.hpp"

#include <vector>
#include <numeric>

using namespace std;

TEST(schedule, solve)
{
  unsigned int stages, nodes;

  stages = 4;
  nodes = 32;

  vector<double> diff{1000.0, 1.0, 1000.0, 1.0};

  auto nps = solve(stages, nodes, diff);

  for(const auto& n : nps)
    printf("%d, ", n);
    
  printf("\nTotal: %d\n", accumulate(nps.begin(), nps.end(), 0));
}

TEST(schedule, to_difficulty)
{
  std::vector<double> perstage_avg = {20.0, 2.0, 10.0, 3.0};
  std::vector<double> expected = {10.0, 1.0, 5.0, 1.5};

  to_difficulty(perstage_avg);

  for (unsigned int i = 0; i < perstage_avg.size(); ++i)
  {
    EXPECT_DOUBLE_EQ(perstage_avg[i], expected[i]);
  }   
  
}

TEST(schedule, assign)
{ 
  int local_stage;

  std::vector<unsigned int> nps = {15, 1, 15, 1};
  std::unordered_map<int, std::unordered_map<int, int>> rank_assignm;
  assign(nps, rank_assignm, &local_stage);

  int stage = 0;
  for (auto& pair : rank_assignm)
  {
    std::unordered_map<int, int> stage_assign = pair.second;
    printf("Stage: %d\nRanks:", stage++);
    for (auto& opair : stage_assign)
    {
      auto rank = opair.second;
      printf(" %d", rank);
    }
    printf("\n");
  }

}