#ifndef EPA_PLACEMENT_SET_H_
#define EPA_PLACEMENT_SET_H_

#include "Placement.hpp"

#include <utility>
#include <string>

class Placement_Set {
public:
  typedef Placement                                        value_type;
  typedef typename std::vector<Placement>::iterator        iterator;
  typedef typename std::vector<Placement>::const_iterator  const_iterator;

  Placement_Set();
  Placement_Set(const unsigned int size);
  Placement_Set(const std::string newick);
  ~Placement_Set();

  // member access
  inline Placement& back() {return placements_.back();};
  inline unsigned int size() {return placements_.size();};
  inline const std::string& newick() const {return newick_;};

  // needs to be in the header
  template<typename ...Args> void emplace_back(Args && ...args)
  {
    placements_.emplace_back(std::forward<Args>(args)...);
  };

  // Iterator Compatability
  inline iterator begin() {return placements_.begin();};
  inline iterator end() {return placements_.end();};
  inline const_iterator begin() const {return placements_.cbegin();};
  inline const_iterator end() const {return placements_.cend();};
  inline const_iterator cbegin() {return placements_.cbegin();};
  inline const_iterator cend() {return placements_.cend();};

private:
  std::vector<Placement> placements_;
  std::string newick_;
};

#endif
