#ifndef EPA_SEQUENCE_H_
#define EPA_SEQUENCE_H_

#include <string>

#include "pllhead.hpp"

class Sequence
{
public:
  Sequence() = default;
  Sequence(std::string header, std::string sequence) : header_(header), sequence_(sequence) {};
  Sequence(const Sequence& s) = default;
  Sequence(Sequence&& s) = default;

  // operator overloads
  Sequence& operator = (const Sequence& s) = default;
  Sequence& operator = (Sequence&& s) = default;

  ~Sequence() = default;

  // member access
  inline const std::string& header() const {return header_;}
  inline const std::string& sequence() const {return sequence_;}

private:
  std::string header_;
  std::string sequence_;

};

#endif
