#pragma once

#include <unordered_map>
#include <vector>

using schedule_type = std::vector< std::vector< int > >;

void to_difficulty( std::vector< double >& perstage_avg );
std::vector< unsigned int > solve( unsigned int stages,
                                   unsigned int nodes,
                                   std::vector< double > const& difficulty_per_stage );
void assign( int const local_rank,
             std::vector< unsigned int >& nodes_per_stage,
             schedule_type& rank_assignm,
             int* local_stage );
void reassign( int const local_rank,
               std::vector< unsigned int >& nodes_per_stage,
               schedule_type& rank_assignm,
               int* local_stage );
