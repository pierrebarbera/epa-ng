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
  ~MSA();
  std::tuple<std::string, std::string> get(const int i) const;

  template <typename T, typename D>
  void append(T header, D sequence);


  iterator begin();
  iterator end();
  const_iterator cbegin();
  const_iterator cend();

  const int num_sites;

private:
  std::vector<Sequence> sequence_list;
};

#endif
