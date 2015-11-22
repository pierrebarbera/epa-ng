#ifndef EPA_PLACEMENT_H_
#define EPA_PLACEMENT_H_

#include <vector>

#include "Sequence.hpp"

class Placement {
public:
  typedef double                                        value_type;
  typedef typename std::vector<double>::iterator        iterator;
  typedef typename std::vector<double>::const_iterator  const_iterator;

  Placement();
  Placement (const unsigned int size, const Sequence & s);
  ~Placement ();
  void set(const unsigned int branch, const double logl);

  // member access
  inline const Sequence& sequence() const {return sequence_;};

  //Iterator Compatability
  iterator begin();
  iterator end();
  const_iterator begin() const;
  const_iterator end() const;
  const_iterator cbegin();
  const_iterator cend();
private:
  /* TODO map might be required, as its likely we have to sort, and just a
    vector would loose the edge locations. plus we may want to immediatly store
    LWR, etc */
  std::vector<double> likelihoods_;
  const Sequence sequence_;
};

#endif
