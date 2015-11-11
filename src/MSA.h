#ifndef EPA_MSA_H_
#define EPA_MSA_H_

#include <vector>

#include "Sequence.h"

class MSA
{
public:
  typedef Sequence                                        value_type;
  typedef typename std::vector<Sequence>::iterator        iterator;
  typedef typename std::vector<Sequence>::const_iterator  const_iterator;

  MSA(const int num_sites);
  MSA();
  // MSA(const MSA& other);
  // MSA operator= (MSA other);
  // MSA operator= (const MSA& other) = delete;
  ~MSA();
  // swap(const MSA& other);
  void append(const std::string& header, const std::string& sequence);
  void append(Sequence s);

  //Member access
  inline int num_sites() const
  {
    return num_sites_;
  }
  Sequence get(const int i) const;

  //Iterator Compatability
  iterator begin();
  iterator end();
  const_iterator cbegin();
  const_iterator cend();



private:
  // Data Members
  int num_sites_;
  std::vector<Sequence> sequence_list_;
};

#endif
