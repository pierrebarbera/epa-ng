#include "Epatest.hpp"

#include "src/schedule.hpp"

#include <vector>
#include <numeric>

using namespace std;

TEST(schedule, solve)
{
  unsigned int stages, nodes;

  stages = 4;
  nodes = 8;

  vector<double> diff{1.0, 2.0, 1.0, 1.0};

  auto nps = solve(stages, nodes, diff);

  for(const auto& n : nps)
    printf("%d, ", n);
    
  printf("\nTotal: %d\n", accumulate(nps.begin(), nps.end(), 0));
}
