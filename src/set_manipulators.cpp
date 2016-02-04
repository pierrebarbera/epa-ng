#include "set_manipulators.hpp"

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

Range superset(Range a, Range b)
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

/* Find duplicate sequences in a MSA and collapse them into one entry that
  holds all respective headers */
void find_collapse_equal_sequences(MSA& msa)
{
  auto end_merge_range = msa.end();
  auto begin_merge_range = end_merge_range;
  for (auto target_iter = msa.begin(); target_iter != end_merge_range; target_iter++)
  {
    auto target = *target_iter;
    begin_merge_range = partition(target_iter + 1, end_merge_range,
    [&target](const Sequence& query)
    {
      return !(target == query);
    });
    // now all sequences in the msa that are equal to the "target" are at the end of the msa
    // whose first element is *begin_merge_range
    for (auto merge_iter = begin_merge_range; merge_iter != end_merge_range; merge_iter++)
      (*target_iter).merge((*merge_iter));

    end_merge_range = begin_merge_range;

  }
  // merging done: all redundant sequences are at the back of msa, starting at end_merge_range
  // cleanup:
  msa.erase(end_merge_range, msa.end());
}
