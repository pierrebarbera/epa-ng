#ifndef EPA_SEQUENCE_H_
#define EPA_SEQUENCE_H_

#include <string>

class Sequence
{
public:
  template<typename T, typename D>
  Sequence(T header, D sequence);
  ~Sequence();
  std::string sequence;
  std::string header;

};

#endif
