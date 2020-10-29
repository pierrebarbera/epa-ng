#pragma once

#include "seq/MSA.hpp"

/* Find duplicate sequences in a MSA and collapse them into one entry that
  holds all respective headers */
void find_collapse_equal_sequences( MSA& msa )
{
  auto end_merge_range   = msa.end();
  auto begin_merge_range = end_merge_range;
  for( auto target_iter = msa.begin(); target_iter != end_merge_range; target_iter++ ) {
    auto target       = *target_iter;
    begin_merge_range = partition( target_iter + 1, end_merge_range,
                                   [&target]( Sequence const& query ) {
                                     return !( target == query );
                                   } );
    // now all sequences in the msa that are equal to the "target" are at the end of the msa
    // whose first element is *begin_merge_range
    for( auto merge_iter = begin_merge_range; merge_iter != end_merge_range; merge_iter++ ) {
      ( *target_iter ).merge( ( *merge_iter ) );
    }

    end_merge_range = begin_merge_range;
  }
  // merging done: all redundant sequences are at the back of msa, starting at end_merge_range
  // cleanup:
  msa.erase( end_merge_range, msa.end() );
}