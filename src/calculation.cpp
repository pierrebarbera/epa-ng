#include "calculation.hpp"

#include <stdexcept>

using namespace std;

void compute_and_set_lwr(PQuery_Set& pqs)
{
  for(auto &pq : pqs)
  {
    auto total = pq.total();
    if (total == 0.0)
      throw runtime_error{"Tried to compute LWR with total = 0"};
    for (auto &p : pq)
      p.lwr(p.likelihood() / total);
  }
}
