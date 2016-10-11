#pragma once

#include <vector>
#include <unordered_map>

void to_difficulty(std::vector<double>& perstage_avg);
std::vector<unsigned int> solve(unsigned int stages, unsigned int nodes, 
            std::vector<double> difficulty_per_stage);
void assign(const int local_rank,
            std::vector<unsigned int>& nodes_per_stage, 
             std::vector<std::vector<int>>& rank_assignm,
            int* local_stage);
void reassign(const int local_rank,
              std::vector<unsigned int>& nodes_per_stage, 
              std::vector<std::vector<int>>& rank_assignm,
              int* local_stage);
