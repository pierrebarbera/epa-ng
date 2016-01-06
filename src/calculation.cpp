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

void discard_by_accumulated_threshold(PQuery_Set& pqs, const double thresh)
{
  // sorting phase
  for (auto &pq : pqs)
  {
    sort(pq.begin(), pq.end(),
      [](const Placement &p_a, const Placement &p_b) -> bool {return p_a.lwr() > p_b.lwr();}
    );
  }
  // accumulation and erasure phase
  for (auto &pq : pqs)
  {
    double sum = 0.0;
    auto pq_iter = pq.begin();
    for (pq_iter = pq.begin(); pq_iter != pq.end(); ++pq_iter)
    {
      sum += pq_iter->lwr();
      if (sum >= thresh)
        break;
    }


    pq.erase(pq_iter, pq.end());
  }
}
