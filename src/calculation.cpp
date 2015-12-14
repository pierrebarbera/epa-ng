#include "calculation.hpp"

#include <stdexcept>
#include <limits>

using namespace std;

void compute_and_set_lwr(PQuery_Set& pqs)
{
  // for(auto &pq : pqs)
  // {
  //   auto total = pq.total();
  //   auto max = pq.max_lh();
  //   total = exp(total - pq.size() * max);
  //
  //   if (total == 0.0)
  //     throw runtime_error{"Tried to compute LWR with total = 0"};
  //   for (auto &p : pq)
  //     p.lwr(exp(p.likelihood() - max) / total);
  // }

  for(auto &pq : pqs)
  {
    double total = 0.0;
    double max=-numeric_limits<double>::infinity();

    for (auto &p : pq)
      if (p.likelihood() > max)
        max = p.likelihood();

    for (auto &p : pq)
      total += exp(p.likelihood() - max);

    for (auto &p : pq)
      p.lwr(exp(p.likelihood() - max) / total);
  }

}
