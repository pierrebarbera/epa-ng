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
  PQuery (const Sequence & s) : sequence_(s) {};
  ~PQuery() = default;

  // needs to be in the header
  template<typename ...Args> void emplace_back(Args && ...args)
  {
    placements_.emplace_back(std::forward<Args>(args)...);
  };

  // member access
  Placement& back() {return placements_.back();};
  const Sequence& sequence() const {return sequence_;};
  unsigned int size() const {return placements_.size();};
  void erase(iterator begin, iterator end) {placements_.erase(begin, end);};

  // Iterator Compatibility
  iterator begin() {return placements_.begin();};
  iterator end() {return placements_.end();};
  const_iterator begin() const {return placements_.cbegin();};
  const_iterator end() const {return placements_.cend();};
  const_iterator cbegin() {return placements_.cbegin();};
  const_iterator cend() {return placements_.cend();};

private:
  std::vector<Placement> placements_;
  const Sequence sequence_;
};

#endif
