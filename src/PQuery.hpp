#ifndef EPA_PQUERY_H_
#define EPA_PQUERY_H_

#include <vector>

#include "Sequence.hpp"
#include "Placement.hpp"

class PQuery {
public:
  typedef Placement                                        value_type;
  typedef typename std::vector<Placement>::iterator        iterator;
  typedef typename std::vector<Placement>::const_iterator  const_iterator;

  PQuery() = delete;
  PQuery (const Sequence & s);
  ~PQuery ();

  // needs to be in the header
  template<typename ...Args> void emplace_back(Args && ...args)
  {
    placements_.emplace_back(std::forward<Args>(args)...);
    // update total
    total_ += placements_.back().likelihood();
  };

  // member access
  inline const Sequence& sequence() const {return sequence_;};
  inline double total() const {return total_;};

  //Iterator Compatability
  iterator begin();
  iterator end();
  const_iterator begin() const;
  const_iterator end() const;
  const_iterator cbegin();
  const_iterator cend();

private:
  std::vector<Placement> placements_;
  const Sequence sequence_;
  double total_;
};

#endif
