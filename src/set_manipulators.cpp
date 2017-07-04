#include "set_manipulators.hpp"

#include <stdexcept>
#include <limits>
#include <algorithm>
#include <iterator>


void split( const Work& src, 
            std::vector<Work>& parts, 
            const unsigned int num_parts)
{
  parts.clear();
  // ensure that there are actually as many parts as specified. We want empty parts to enable null messages
  parts.resize(num_parts);

  const size_t ext_size = (src.size() - (src.size() % num_parts)) + num_parts;
  const size_t chunk_size = ext_size / num_parts;

  size_t i = 0;
  size_t bucket = 0;
  for (auto it : src) {
    parts[bucket].add(it);
    if(++i % chunk_size == 0) bucket++;
  }
}


void merge(Work& dest, const Work& src)
{
  auto prev_branch_id = (*src.begin()).branch_id + 1;
  for (auto it : src)  {
    const auto branch_id = it.branch_id;
    if (prev_branch_id != branch_id) {
      dest[branch_id];
      dest[branch_id].insert( dest.at(branch_id).end(),
                              src.at(branch_id).begin(), 
                              src.at(branch_id).end()
                            );
    }
    prev_branch_id = branch_id;
  }
}

void merge(Timer<>& dest, const Timer<>& src)
{
  dest.insert(dest.end(), src.begin(), src.end());
}

void compute_and_set_lwr(Sample<Placement>& sample)
{
  for (auto &pq : sample) {
    double total = 0.0;
    double max=-std::numeric_limits<double>::infinity();

    double entropy = 0.0;

    // find the maximum
    for (auto &p : pq) {
      if (p.likelihood() > max) {
        max = p.likelihood();
      }
    }

    // sum up the distances to the max
    for (auto &p : pq) {
      total += exp(p.likelihood() - max);
    }

    // normalize the distances
    for (auto &p : pq) {
      double lwr = exp(p.likelihood() - max) / total;
      p.lwr(lwr);
      // compute the shannon entropy of the query (in nats)
      entropy -= (lwr * log(lwr));
    }

    pq.entropy(entropy);
  }
}

static void sort_by_lwr(PQuery<Placement>& pq)
{
  sort(pq.begin(), pq.end(),
    [](const Placement &p_a, const Placement &p_b) -> bool {
      return p_a.lwr() > p_b.lwr();
    }
  );
}

void discard_bottom_x_percent(Sample<Placement>& sample, const double x)
{
  if (x < 0.0 || x > 1.0) {
    throw std::range_error{"x is not a percentage (outside of [0,1])"};
  }

  for (auto &pq : sample) {
    auto num_keep = (int)ceil((1.0 - x) * static_cast<double>(pq.size()));
    sort_by_lwr(pq);
    auto erase_iter = pq.begin();
    advance(erase_iter, num_keep);
    pq.erase(erase_iter, pq.end());
  }
}

void discard_by_support_threshold(Sample<Placement>& sample, 
                                  const double thresh, 
                                  const unsigned int min, 
                                  const unsigned int max)
{
  if (thresh < 0.0 || thresh > 1.0){
    throw std::range_error{"thresh is not a valid likelihood weight ratio (outside of [0,1])"};
  }

  if (min < 1) {
    throw std::range_error{"Filter min cannot be smaller than 1!"};
  }

  if (min > max) {
    throw std::range_error{"Filter min cannot be smaller than max!"};
  }

  // static_assert(std::is_array<PQuery>::value,
  //                 "PQuery not array type");

  for (auto &pq : sample) {
    auto erase_iter = partition(
      pq.begin(), 
      pq.end(),
      [thresh](Placement &p) -> bool {
        return (p.lwr() > thresh);
      }
    );

    if ( distance(erase_iter, (pq.begin() + min - 1)) > 0 ) {
      erase_iter = pq.begin() + min - 1;
    }

    if ( pq.size() <= max && distance((pq.begin() + max - 1), erase_iter) > 0 ) {
      erase_iter = pq.begin() + max - 1;
    }

    pq.erase(erase_iter, pq.end());
  }
}

void discard_by_accumulated_threshold(Sample<Placement>& sample, 
                                      const double thresh,
                                      const unsigned int min, 
                                      const unsigned int max)
{
  if (thresh < 0.0 || thresh > 1.0) {
    throw std::range_error{"thresh is not a valid likelihood weight ratio (outside of [0,1])"};
  }

  if (min < 1) {
    throw std::range_error{"Filter min cannot be smaller than 1!"};
  }
  
  if (min > max) {
    throw std::range_error{"Filter min cannot be smaller than max!"};
  }

  // sorting phase
  for (auto &pq : sample) {
    sort_by_lwr(pq);
  }

  // accumulation and erasure phase
  for (auto &pq : sample) {
    double sum = 0.0;

    auto pq_iter = pq.begin();
    const auto max_iter = ( distance(pq_iter + max - 1, pq.end()) > 0 ) ? pq_iter + max - 1 : pq.end();
    for (pq_iter = pq.begin(); pq_iter != max_iter && sum < thresh; ++pq_iter) {
      // sum up until threshold is passed. if we abort before it is passed, we would have the possibility of
      // empty lists
      sum += pq_iter->lwr();
    }

    auto to_add = distance(pq_iter, pq.begin() + min - 1);
    if (to_add > 0) {
      advance(pq_iter, to_add);
    }

    pq.erase(pq_iter, pq.end());
  }
}

Range superset(Range a, Range b)
{
  if (b.begin < a.begin) {
    a.begin = b.begin;
  }
  if (b.span > a.span) {
    a.span = b.span;
  }
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
Range get_valid_range(std::string sequence)
{
  unsigned int lower = 0;
  unsigned int upper = sequence.length();

  while(sequence.c_str()[lower] == '-') {
    lower++;
  }

  while(sequence.c_str()[upper - 1] == '-') {
    upper--;
  }

  return Range(lower, upper - lower);
}

/* Find duplicate sequences in a MSA and collapse them into one entry that
  holds all respective headers */
void find_collapse_equal_sequences(MSA& msa)
{
  auto end_merge_range = msa.end();
  auto begin_merge_range = end_merge_range;
  for (auto target_iter = msa.begin(); target_iter != end_merge_range; target_iter++) {
    auto target = *target_iter;
    begin_merge_range = partition(target_iter + 1, end_merge_range,
    [&target](const Sequence& query)
    {
      return !(target == query);
    });
    // now all sequences in the msa that are equal to the "target" are at the end of the msa
    // whose first element is *begin_merge_range
    for (auto merge_iter = begin_merge_range; merge_iter != end_merge_range; merge_iter++) {
      (*target_iter).merge((*merge_iter));
    }

    end_merge_range = begin_merge_range;

  }
  // merging done: all redundant sequences are at the back of msa, starting at end_merge_range
  // cleanup:
  msa.erase(end_merge_range, msa.end());
}
