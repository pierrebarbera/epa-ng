#ifndef EPA_CALCULATION_H_
#define EPA_CALCULATION_H_

#include "PQuery_Set.hpp"

void compute_and_set_lwr(PQuery_Set& pqs);
void discard_by_support_threshold(PQuery_Set& pqs, const double thresh);


#endif
