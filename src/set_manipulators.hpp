#pragma once

#include <vector>
#include <unordered_map>
#include <limits>

#include "sample/Sample.hpp"
#include "util/Timer.hpp"
#include "util/Options.hpp"
#include "util/logging.hpp"
#include "seq/MSA.hpp"
#include "core/Work.hpp"

#include "sample/Device_Sample.hpp"
#include <thrust/device_vector.h>
#include <thrust/device_ptr.h>

using pq_citer_t  = PQuery::const_iterator;
using pq_iter_t   = PQuery::iterator;

void merge(Work& dest, const Work& src);
void merge(Timer<>& dest, const Timer<>& src);

void sort_by_lwr(PQuery& pq);
void sort_by_logl(PQuery& pq);
void compute_and_set_lwr( Sample& sample );
pq_iter_t until_top_percent( PQuery& pq,
                              const double x);
void discard_bottom_x_percent(Sample& sample, const double x);
void discard_by_support_threshold(Sample& sample,
                                  const double thresh,
                                  const size_t min=1,
                                  const size_t max=std::numeric_limits<size_t>::max());


pq_iter_t until_accumulated_reached(  PQuery& pq,
                                      const double thresh,
                                      const size_t min,
                                      const size_t max);

pq_iter_t until_accumulated_reached(  PQuery& pq,
                                      const double thresh);

void discard_by_accumulated_threshold(Sample& sample,
                                      const double thresh,
                                      const size_t min=1,
                                      const size_t max=std::numeric_limits<size_t>::max());
void filter(Sample& sample, const Options& options);
void find_collapse_equal_sequences(MSA& msa);

/**
 * collapses PQuerys with the same ID inside a Sample into one
 */
void collapse(Sample& sample);

/**
 * dummy collapse function for default case of objects where
 * collapsing doesnt make sense
 */
template<class T>
void collapse(T&)
{ }

/**
 * special split function that Splits samples in buckets according to the global sequence ID
 * of their PQueries. The goal is to have them split such that each aggregate node gets their
 * correct set of sequence results (even if that part is empty, which constitutes a null-message)
 *
 * @param
 * @param
 * @param
 */
// template<class T>
void split( const Sample& src,
            std::vector<Sample>& parts,
            const unsigned int num_parts);

/**
 * Splits a <src> into <num_parts> number of equally sized <parts>.
 * <parts> is expected to be empty.
 *
 * REQUIRES T TO BE:
 *   Iterable
 *   has size()
 *   has insert()
 *
 * @param src    Container to split
 * @param parts     resulting parts vector
 * @param num_parts number of parts
 */
template <class T>
void split( const T& src,
            std::vector<T>& parts,
            const unsigned int num_parts)
{
  parts.clear();
  unsigned int chunk_size = ceil(src.size() / static_cast<double>(num_parts));
  auto move_begin = src.begin();

  while (move_begin < src.end()) {
    auto move_end = std::min(move_begin + chunk_size, src.end());
    parts.emplace_back();
    parts.back().insert(move_begin, move_end);
    move_begin = move_end;
  }
}

void split( const Work& source,
            std::vector<Work>& parts,
            const unsigned int num_parts);

/**
  Merges a Sample <src> into a Sample <dest>. Leaves <src> intact.
*/
// template<class T>
void merge(Sample& dest, const Sample& src);

/**
  Merges a Sample <src> into a Sample <dest>. src here is an rvalue,
  and thus the elements are moved instead of copied
*/
// template<class T>
void merge(Sample& dest, Sample&& src);

template<class T>
void merge(T& dest, std::vector<T>&& parts)
{
  for (auto& p : parts) {
    merge(dest, std::move(p));
  }
}

template <class T>
void merge(T& dest, std::vector<T>& parts)
{
  for (const auto& p : parts) {
    merge(dest, p);
  }
}

template <class T>
void merge(std::vector<T>& dest, const std::vector<T>& parts)
{
  dest.insert( std::end(dest), std::begin(parts), std::end(parts) );
}
