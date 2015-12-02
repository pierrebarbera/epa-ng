#ifndef EPA_PQUERY_SET_H_
#define EPA_PQUERY_SET_H_

#include "PQuery.hpp"

#include <utility>
#include <string>

class PQuery_Set {
public:
  typedef PQuery                                        value_type;
  typedef typename std::vector<PQuery>::iterator        iterator;
  typedef typename std::vector<PQuery>::const_iterator  const_iterator;

  PQuery_Set() = default;
  PQuery_Set(const std::string newick);
  ~PQuery_Set() = default;

  // member access
  inline PQuery& back() {return pquerys_.back();};
  inline unsigned int size() {return pquerys_.size();};
  inline const std::string& newick() const {return newick_;};

  // needs to be in the header
  template<typename ...Args> void emplace_back(Args && ...args)
  {
    pquerys_.emplace_back(std::forward<Args>(args)...);
  };

  // Iterator Compatability
  inline iterator begin() {return pquerys_.begin();};
  inline iterator end() {return pquerys_.end();};
  inline const_iterator begin() const {return pquerys_.cbegin();};
  inline const_iterator end() const {return pquerys_.cend();};
  inline const_iterator cbegin() {return pquerys_.cbegin();};
  inline const_iterator cend() {return pquerys_.cend();};

private:
  std::vector<PQuery> pquerys_;
  std::string newick_;
};

#endif
