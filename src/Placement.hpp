#ifndef EPA_PLACEMENT_H_
#define EPA_PLACEMENT_H_

#include <vector>

#include "Sequence.hpp"

class Placement {
public:
  Placement (const unsigned int size, const Sequence * s);
  ~Placement ();
  void set(const unsigned int branch, const double logl);
private:
  std::vector<double> likelihoods_;
  const Sequence * sequence;
};

#endif
