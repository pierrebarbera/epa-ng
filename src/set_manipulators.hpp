#pragma once

#include "PQuery_Set.hpp"
#include "Range.hpp"

void compute_and_set_lwr(PQuery_Set& pqs);
void discard_by_support_threshold(PQuery_Set& pqs, const double thresh);
void discard_by_accumulated_threshold(PQuery_Set& pqs, const double thresh);
Range superset(Range a, Range b);
Range get_valid_range(std::string sequence);
