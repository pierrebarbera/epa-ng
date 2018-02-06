#pragma once

#include <algorithm>

#include "sample/Sample.hpp"
#include "util/Options.hpp"
#include "set_manipulators.hpp"

inline void dynamic_heuristic(Sample<Placement>& sample,
                              const Options& options)
{
  compute_and_set_lwr(sample);
  discard_by_accumulated_threshold( sample, 
                                    options.prescoring_threshold,
                                    options.filter_min,
                                    options.filter_max);
}

inline void fixed_heuristic(Sample<Placement>& sample,
                            const Options& options)
{
  compute_and_set_lwr(sample);
  discard_bottom_x_percent(sample,
                          (1.0 - options.prescoring_threshold));
}

static inline void sort_by_logl(PQuery<Placement>& pq)
{
  std::sort(pq.begin(), pq.end(),
    [](const Placement &lhs, const Placement &rhs) -> bool {
      return lhs.likelihood() > rhs.likelihood();
    }
  );
}

inline void baseball_heuristic( Sample<Placement>& sample,
                                const Options& options)
{
  (void) options;
  // strike_box: logl delta, keep placements within this many logl units from the best
  const double strike_box = 3;
  // max_strikes: number of additional branches to add after strike box is full
  const size_t max_strikes = 6;
  // max_pitches: absolute maximum of candidates to select
  const size_t max_pitches = 40;
  
  #ifdef __OMP
  #pragma omp parallel for schedule(dynamic)
  #endif
  for (size_t i = 0; i < sample.size(); ++i) {
    auto &pq = sample[i];
    assert(pq.size());
    // sort placements by likelihood (descending)
    sort_by_logl(pq);
    // keep any placements that are within strike box of the best
    const double best_logl = pq[0].likelihood();
    const double thresh = best_logl - strike_box;
    // get first element not within strike box
    auto keep_iter = std::find_if(pq.begin(), pq.end(),
      [thresh](const auto& p){
        return (p.likelihood() < thresh);
      }
    );

    const auto hits = std::distance(pq.begin(), keep_iter);

    // ensure we keep no more than max_pitches
    size_t to_add = std::min(max_pitches - hits, max_strikes);

    std::advance(keep_iter, to_add);

    // erase any other placements that were not selected
    pq.erase(keep_iter, pq.end());
  }
}

inline void apply_heuristic(Sample<Placement>& sample,
                            const Options& options)
{
  if (options.baseball) {
    baseball_heuristic(sample, options);
  } else if (options.prescoring_by_percentage) {
    dynamic_heuristic(sample, options);
  } else {
    fixed_heuristic(sample, options);    
  }
}
