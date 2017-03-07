#include "set_manipulators.hpp"

#include <stdexcept>
#include <limits>
#include <algorithm>
#include <iterator>

using namespace std;

/**
  Splits a Sample <source> into <parts> according to <split_map>.

  <split_map> is structured as such: splitmap[part_id] = {<set of sequence id's that will be moved to the part>}
  assumes <parts> to be empty / that its content won't be missed

*/
void split(Sample& src, vector<Sample>& parts, const vector<vector<unsigned int>>& split_map)
{
  parts.clear();
  for (const auto& part_move_list : split_map)
  {
    parts.push_back(Sample());
    // move all instances of pquery specified in split_map for the current part from src to parts[part_id]
    for (auto sequence_to_move : part_move_list)
    {
      // using already implemented == behaviour of PQuery to find the correct pquery to move
      PQuery query(sequence_to_move);
      auto to_move = find(src.begin(), src.end(), query);

      if (to_move != src.end())
      {
        // auto move_to = parts.back().end();
        // move(to_move, to_move + 1, move_to);
        // move the entry (leaves original entry in undefined but valid state)
        parts.back().push_back(move(*to_move));

        // remove the original entry
        src.erase(to_move, to_move + 1);
      }
      // TODO fail if sequence not found?
    }
  }
  assert(parts.size() == split_map.size());
}

/**
 * special split function that Splits samples in buckets according to the global sequence ID
 * of their PQueries. The goal is to have them split such that each aggregate node gets their
 * correct set of sequence results (even if that part is empty, which constitutes a null-message)
 * 
 * @param
 * @param
 * @param
 */
void split(const Sample& src, std::vector<Sample>& parts, const unsigned int num_parts)
{
  parts.clear();
  // ensure that there are actually as many parts as specified. We want empty parts to enable null messages
  parts.resize(num_parts);

  for (auto& pq : src)
  {
    const auto bucket = pq.sequence_id() % num_parts;
    parts[bucket].push_back(pq);
  }

}

void split(const Work& src, std::vector<Work>& parts, const unsigned int num_parts)
{
  parts.clear();
  // ensure that there are actually as many parts as specified. We want empty parts to enable null messages
  parts.resize(num_parts);

  const size_t ext_size = (src.size() - (src.size() % num_parts)) + num_parts;
  const size_t chunk_size = ext_size / num_parts;

  size_t i = 0;
  size_t bucket = 0;
  for (auto it : src)
  {
    parts[bucket].add(it);
    if(++i % chunk_size == 0) bucket++;
  }
}

/**
  Merges a Sample <src> into a Sample <dest>. Leaves <src> intact.
*/
void merge(Sample& dest, const Sample& src)
{
  // merge in every source pquery...
  for (const auto& pquery : src)
  {
    // ... by checking if its sequence already exists in destination
    auto input_iter = find(dest.begin(), dest.end(), pquery);
    // if not, create a record
    if (input_iter == dest.end())
    {
      dest.emplace_back(pquery.sequence_id());
      input_iter = --(dest.end());
    }
    // then concat their vectors
    input_iter->insert(input_iter->end(), pquery.begin(), pquery.end());
  }
}

void merge(Work& dest, const Work& src)
{
  auto prev_branch_id = (*src.begin()).branch_id + 1;
  for (auto it : src) 
  {
    const auto branch_id = it.branch_id;
    if (prev_branch_id != branch_id)
    {
      dest[branch_id];
      dest[branch_id].insert( dest.at(branch_id).end(),
                              src.at(branch_id).begin(), 
                              src.at(branch_id).end()
                            );
    }
    prev_branch_id = branch_id;
  }
}

void merge(Timer& dest, const Timer& src)
{
  dest.insert(dest.end(), src.begin(), src.end());
}

void compute_and_set_lwr(Sample& sample)
{
  for (auto &pq : sample)
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

static void sort_by_lwr(PQuery& pq)
{
  sort(pq.begin(), pq.end(),
    [](const Placement &p_a, const Placement &p_b) -> bool {return p_a.lwr() > p_b.lwr();}
  );
}

void discard_bottom_x_percent(Sample& sample, const double x)
{
  if (x < 0.0 || x > 1.0)
    throw range_error{"x is not a percentage (outside of [0,1])"};
  for (auto &pq : sample)
  {
    auto num_keep = (int)ceil((1.0 - x) * (double)pq.size());
    sort_by_lwr(pq);
    auto erase_iter = pq.begin();
    advance(erase_iter, num_keep);
    pq.erase(erase_iter, pq.end());
  }
}

void discard_by_support_threshold(Sample& sample, const double thresh, 
                                  const unsigned int min, const unsigned int max)
{
  if (thresh < 0.0 || thresh > 1.0)
    throw range_error{"thresh is not a valid likelihood weight ratio (outside of [0,1])"};

  if (min < 1)
    throw range_error{"Filter min cannot be smaller than 1!"};

  if (min > max)
    throw range_error{"Filter min cannot be smaller than max!"};

  // static_assert(std::is_array<PQuery>::value,
  //                 "PQuery not array type");

  for (auto &pq : sample)
  {
    auto erase_iter = partition(
      pq.begin(), 
      pq.end(),
      [thresh](Placement &p) -> bool {
        return (p.lwr() > thresh);
      }
    );

    if ( distance(erase_iter, (pq.begin() + min - 1)) > 0 )
      erase_iter = pq.begin() + min - 1;

    if ( pq.size() <= max && distance((pq.begin() + max - 1), erase_iter) > 0 )
      erase_iter = pq.begin() + max - 1;

    pq.erase(erase_iter, pq.end());
  }
}

void discard_by_accumulated_threshold(Sample& sample, const double thresh,
                                  const unsigned int min, const unsigned int max)
{
  if (thresh < 0.0 || thresh > 1.0)
    throw range_error{"thresh is not a valid likelihood weight ratio (outside of [0,1])"};

  if (min < 1)
    throw range_error{"Filter min cannot be smaller than 1!"};
  
  if (min > max)
    throw range_error{"Filter min cannot be smaller than max!"};

  // sorting phase
  for (auto &pq : sample)
    sort_by_lwr(pq);

  // accumulation and erasure phase
  for (auto &pq : sample)
  {
    double sum = 0.0;


    auto pq_iter = pq.begin();
    const auto max_iter = ( distance(pq_iter + max - 1, pq.end()) > 0 ) ? pq_iter + max - 1 : pq.end();
    for (pq_iter = pq.begin(); pq_iter != max_iter && sum < thresh; ++pq_iter)
    {
      // sum up until threshold is passed. if we abort before it is passed, we would have the possibility of
      // empty lists
      sum += pq_iter->lwr();
      
    }

    auto to_add = distance(pq_iter, pq.begin() + min - 1);
    if (to_add > 0)
      advance(pq_iter, to_add);

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
