#include "set_manipulators.hpp"

#include <stdexcept>
#include <limits>
#include <algorithm>
#include <iterator>
#include <cmath>

void collapse(Sample& sample)
{
  const auto invalid = std::numeric_limits<
    typename PQuery::seqid_type>::max();

  std::unordered_map< size_t, std::vector<size_t> > collapse_set;

  // build map of all pqueries
  for (size_t i = 0; i < sample.size(); ++i) {
    const auto& pq = sample[i];
    collapse_set[pq.sequence_id()].emplace_back(i);
  }

  // find all cases of puplicate entries and merge hem
  for (auto& pair : collapse_set) {
    auto pqlist = pair.second;
    // duplicate!
    if (pqlist.size() > 1) {
      // move entries from duplicate to original
      auto& dest = sample[pqlist[0]].data();
      for (size_t i = 1; i < pqlist.size(); ++i) {
        auto& src = sample[pqlist[i]].data();
        dest.reserve(dest.size() + src.size());
        std::move(std::begin(src), std::end(src), std::back_inserter(dest));
        // mark invalid in original sample
        sample[pqlist[i]].sequence_id(invalid);
      }
    }
  }

  // clear the original sample of invalid pqueries
  sample.erase(
    std::remove_if( std::begin(sample),
                    std::end(sample),
                    [invalid = invalid](auto& e){
                      return e.sequence_id() == invalid;}),
    std::end(sample) );

}

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

void split( const Sample& src,
            std::vector<Sample>& parts,
            const unsigned int num_parts)
{
  parts.clear();
  // ensure that there are actually as many parts as specified. We want empty parts to enable null messages
  parts.resize(num_parts);

  for (auto& pq : src) {
    const auto bucket = pq.sequence_id() % num_parts;
    parts[bucket].push_back(pq);
  }

}

void merge(Work& dest, const Work& src)
{
  if ( not src.empty() ) {
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
}

void merge(Sample& dest, const Sample& src)
{
  // merge in every source pquery...
  for (const auto& pquery : src) {
    // ... by checking if its sequence already exists in destination
    auto input_iter = find(dest.begin(), dest.end(), pquery);
    // if not, create a record
    if (input_iter == dest.end()) {
      dest.emplace_back(pquery.sequence_id(), pquery.header());
      input_iter = --(dest.end());
    }
    // then concat their vectors
    input_iter->insert(input_iter->end(), pquery.begin(), pquery.end());
  }
}

/**
  Merges a Sample <src> into a Sample <dest>. src here is an rvalue,
  and thus the elements are moved instead of copied
*/
// template<class T>
void merge(Sample& dest, Sample&& src)
{
  for (auto& pquery : src) {
    // create new record
    dest.emplace_back(std::move(pquery));
  }
}

// void merge(Sample& dest, Device_Sample&& src)
// {
//   thrust::host_vector<Device_PQuery> tmp = src.pquerys_;
//   for (auto& pquery : tmp) {
//     // create new record
//     dest.emplace_back(std::move(pquery));
//   }
// }

void merge(Timer<>& dest, const Timer<>& src)
{
  dest.insert(dest.end(), src.begin(), src.end());
}

void compute_and_set_lwr(Sample& sample)
{
  #ifdef __OMP
  #pragma omp parallel for schedule(dynamic)
  #endif
  for (size_t j = 0; j < sample.size(); ++j) {
    auto &pq = sample[j];
    double total = 0.0;

    // find the maximum
    auto max = std::max_element(pq.begin(), pq.end(),
      [](const Placement& lhs, const Placement& rhs){
        return (lhs.likelihood() < rhs.likelihood());
      }
    )->likelihood();

    std::vector<double> exp_lh(pq.size());
    //get the distances to the max
    for (size_t i = 0; i < pq.size(); ++i) {
      exp_lh[i] = std::exp(pq[i].likelihood() - max);
      total += exp_lh[i];
    }

    // normalize the distances
    for (size_t i = 0; i < pq.size(); ++i) {
      pq[i].lwr(exp_lh[i] / total);
    }
  }
}

void sort_by_lwr(PQuery& pq)
{
  sort(pq.begin(), pq.end(),
    [](const Placement &p_a, const Placement &p_b) -> bool {
      return p_a.lwr() > p_b.lwr();
    }
  );
}

void sort_by_logl(PQuery& pq)
{
  std::sort(pq.begin(), pq.end(),
    [](const Placement &lhs, const Placement &rhs) -> bool {
      return lhs.likelihood() > rhs.likelihood();
    }
  );
}

pq_iter_t until_top_percent( PQuery& pq,
                              const double x)
{
  sort_by_lwr(pq);
  auto num_keep = static_cast<size_t>(ceil(x * static_cast<double>(pq.size())));
  auto iter = pq.begin();
  advance(iter, num_keep);
  return iter;
}

pq_iter_t until_accumulated_reached( PQuery& pq,
                                      const double thresh,
                                      const size_t min,
                                      const size_t max)
{
  sort_by_lwr(pq);

  double sum = 0.0;
  size_t num_summed = 0;

  auto pq_iter = pq.begin();
  for (pq_iter = pq.begin(); num_summed < max and sum < thresh; ++pq_iter, ++num_summed) {
    // sum up until threshold is passed. if we abort before it is passed, we would have the possibility of
    // empty lists
    sum += pq_iter->lwr();
  }

  auto to_add = distance(pq_iter, pq.begin() + min - 1);
  if (to_add > 0) {
    advance(pq_iter, to_add);
  }

  return pq_iter;
}

pq_iter_t until_accumulated_reached(  PQuery& pq,
                                      const double thresh)
{
  return until_accumulated_reached( pq, thresh, 1, std::numeric_limits<size_t>::max() );
}

void discard_bottom_x_percent(Sample& sample, const double x)
{
  if (x < 0.0 || x > 1.0) {
    throw std::range_error{"x is not a percentage (outside of [0,1])"};
  }

  #ifdef __OMP
  #pragma omp parallel for schedule(dynamic)
  #endif
  for (size_t i = 0; i < sample.size(); ++i) {
    auto &pq = sample[i];
    auto erase_iter = until_top_percent(pq, 1.0 - x);
    pq.erase(erase_iter, pq.end());
  }
}

void discard_by_support_threshold(Sample& sample,
                                  const double thresh,
                                  const size_t min,
                                  const size_t max)
{
  if (thresh < 0.0 or thresh > 1.0){
    throw std::range_error{"thresh is not a valid likelihood weight ratio (outside of [0,1])"};
  }
  if (min < 1) {
    throw std::range_error{"Filter min cannot be smaller than 1!"};
  }

  #ifdef __OMP
  #pragma omp parallel for schedule(dynamic)
  #endif
  for (size_t i = 0; i < sample.size(); ++i) {
    auto &pq = sample[i];
    sort_by_lwr(pq);
    auto erase_iter = partition(
      pq.begin(),
      pq.end(),
      [thresh](Placement &p) -> bool {
        return (p.lwr() > thresh);
      }
    );

    const auto num_kept = static_cast<size_t>(distance(pq.begin(), erase_iter));

    if ( num_kept < min ) {
      const auto to_add = min - num_kept;
      std::advance(erase_iter, to_add);
    }

    if ( max and num_kept > max ) {
      const auto to_remove = num_kept - max;
      std::advance(erase_iter, -to_remove);
    }

    pq.erase(erase_iter, pq.end());
  }
}

void discard_by_accumulated_threshold(Sample& sample,
                                      const double thresh,
                                      const size_t min,
                                      const size_t max)
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

  // accumulation and erasure phase
  #ifdef __OMP
  #pragma omp parallel for schedule(dynamic)
  #endif
  for (size_t i = 0; i < sample.size(); ++i) {
    auto &pq = sample[i];

    auto pq_iter = until_accumulated_reached(pq, thresh, min, max);

    pq.erase(pq_iter, pq.end());
  }
}

void filter(Sample& sample, const Options& options)
{
    if (options.acc_threshold) {
      LOG_DBG << "Filtering output by accumulated threshold: " << options.support_threshold << std::endl;
      discard_by_accumulated_threshold( sample,
                                        options.support_threshold,
                                        options.filter_min,
                                        options.filter_max);
    } else {
      LOG_DBG << "Filtering output placements below threshold: " << options.support_threshold << std::endl;
      discard_by_support_threshold( sample,
                                    options.support_threshold,
                                    options.filter_min,
                                    options.filter_max);
    }
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
