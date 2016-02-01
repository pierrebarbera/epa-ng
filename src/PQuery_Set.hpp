#pragma once

#include "PQuery.hpp"

#include <utility>
#include <string>

class PQuery_Set {
public:
  typedef PQuery                                        value_type;
  typedef typename std::vector<PQuery>::iterator        iterator;
  typedef typename std::vector<PQuery>::const_iterator  const_iterator;

  PQuery_Set() = default;
  PQuery_Set(const std::string newick) : newick_(newick) {};
  ~PQuery_Set() = default;

  // member access
  PQuery& back() {return pquerys_.back();};
  unsigned int size() {return pquerys_.size();};
  const std::string& newick() const {return newick_;};

  // needs to be in the header
  template<typename ...Args> void emplace_back(Args && ...args)
  {
    pquerys_.emplace_back(std::forward<Args>(args)...);
  };

  // Iterator Compatibility
  iterator begin() {return pquerys_.begin();};
  iterator end() {return pquerys_.end();};
  const_iterator begin() const {return pquerys_.cbegin();};
  const_iterator end() const {return pquerys_.cend();};
  const_iterator cbegin() {return pquerys_.cbegin();};
  const_iterator cend() {return pquerys_.cend();};

private:
  // TODO if this is indexed by unique query indexes it could be good for
  // resukt retrieval in mp*
  std::vector<PQuery> pquerys_;
  std::string newick_;
};
