#ifndef EPA_PQUERY_H_
#define EPA_PQUERY_H_

#include <vector>

#include "Sequence.hpp"

class PQuery {
public:
  typedef double                                        value_type;
  typedef typename std::vector<double>::iterator        iterator;
  typedef typename std::vector<double>::const_iterator  const_iterator;

  PQuery();
  PQuery (const unsigned int size, const Sequence & s);
  ~PQuery ();
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
  /* TODO change to vector of pquery class, to make it sortable, */
  std::vector<double> likelihoods_;
  const Sequence sequence_;
};

#endif
