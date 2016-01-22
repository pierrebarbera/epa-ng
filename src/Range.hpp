#pragma once

class Range {
public:
  Range (unsigned int begin, unsigned int span)
    : begin(begin), span(span) {};
  Range() = default;
  ~Range () = default;

  unsigned int begin;
  unsigned int span;

private:

};
