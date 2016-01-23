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

    // find the maximum
    for (auto &p : pq)
      if (p.likelihood() > max)
        max = p.likelihood();

    // sum up the distances to the max
    for (auto &p : pq)
      total += exp(p.likelihood() - max);

    // normalize the distances
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
    for (pq_iter = pq.begin(); pq_iter != pq.end() && sum < thresh; ++pq_iter)
      sum += pq_iter->lwr();
      // sum up until threshold is passed. if we abort before it is passed, we would have the possibility of
      // empty lists

    pq.erase(pq_iter, pq.end());
  }
}

inline Range superset(Range a, Range b)
{
  if (b.begin < a.begin)
    a.begin = b.begin;
  if (b.span > a.span)
    a.span = b.span;
  return a;
}

/*  Returns the range of a sequence outside of which there are ONLY indel characters.
 *  Range starts at the first valid position and ends after <span> characters, where
 *  begin + span is the first element not included in the valid range.
 *  Example:
 *  -  -  -  A  T  A  G  C  T  -  -
 *  0  1  2  3  4  5  6  7  8  9 10
 *  Output: (3,6)
 */
Range get_valid_range(string sequence)
{
  unsigned int lower = 0;
  unsigned int upper = sequence.length();

  while(sequence.c_str()[lower] == '-')
    lower++;

  while(sequence.c_str()[upper - 1] == '-')
    upper--;

  return Range(lower, upper - lower);
}
