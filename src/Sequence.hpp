#ifndef EPA_SEQUENCE_H_
#define EPA_SEQUENCE_H_

#include <string>

class Sequence
{
public:
  Sequence(const Sequence& s);
  Sequence(std::string header, std::string sequence);
  ~Sequence();

  std::string header() const
  {
    return header_;
  }

  const std::string& sequence() const
  {
    return sequence_;
  }

private:
  std::string header_;
  std::string sequence_;

};

#endif
