#include "calculation.hpp"

#include <stdexcept>
#include <limits>
#include <algorithm>

using namespace std;

void compute_and_set_lwr(PQuery_Set& pqs)
{
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

void discard_by_support_threshold(PQuery_Set& pqs, const double thresh)
{
  for (auto &pq : pqs)
  {
    auto erase_iter = partition(pq.begin(), pq.end(),
    [thresh](Placement &p) -> bool {return (p.lwr() > thresh);});
    pq.erase(erase_iter, pq.end());
  }
}
