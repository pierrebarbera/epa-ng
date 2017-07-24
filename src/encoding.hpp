#pragma once

#include <algorithm>
#include <cassert>
#include <cmath>

#include "Matrix.hpp"
#include "maps.hpp"

class FourBit
{
  // shorthand
  using uchar = unsigned char;
private:
  static constexpr uchar NONE_CHAR = 'Z';

  uchar pack_(const uchar lhs, const uchar rhs)
  {
    assert(lhs < 16 and rhs < 16);
    return (lhs << 4) | rhs;
  }

  std::pair<uchar, uchar> unpack_(const uchar c)
  {
    const uchar upper_mask = 0b11110000;
    const uchar lower_mask = 0b00001111;
    return {(c & upper_mask) >> 4, c & lower_mask};
  }

public:
  FourBit()
    : to_fourbit_(128, 128, 16)
  {
    static_assert(NT_MAP_SIZE == 16, "Weird NT map size, go adjust encoder code!");

    // both chars are valid characters:
    const auto map_size = static_cast<uchar>(NT_MAP_SIZE);
    for (uchar i = 0; i < map_size; ++i) {
      assert(NT_MAP[i] == std::toupper(NT_MAP[i]));
      for (uchar j = 0; j < map_size; ++j) {
        uchar packed_char = pack_(i,j);

        // uppercase
        auto row = NT_MAP[i];
        auto col = NT_MAP[j];
        to_fourbit_.at(row, col) = packed_char;
        // lowercase
        row = std::tolower(NT_MAP[i]);
        to_fourbit_.at(row, col) = packed_char;

        col = std::tolower(NT_MAP[j]);
        to_fourbit_.at(row, col) = packed_char;

        row = NT_MAP[i];
        to_fourbit_.at(row, col) = packed_char;
      }
    }

    // valid + padding
    for (uchar i = 0; i < map_size; ++i) {
      uchar packed_char = pack_(i,0);

      // uppercase
      auto row = NT_MAP[i];
      auto col = NONE_CHAR;
      to_fourbit_.at(row, col) = packed_char;

      row = std::tolower(NT_MAP[i]);
      to_fourbit_.at(row, col) = packed_char;
    }
      
  }
  ~FourBit() = default;

  inline size_t packed_size(const size_t size)
  {
    return std::ceil(size / 2.0);
  }

  // conversion functions
  inline std::basic_string<char> to_fourbit(const std::string& s)
  {
    const size_t p_size = packed_size(s.size());
    std::basic_string<char> res;
    res.reserve(p_size);

    size_t i = 0;
    for (; i + 1 < s.size(); i += 2) {
      res.push_back(to_fourbit_.at(s[i], s[i+1u]));
    }

    // original string size not divisible by 2: trailing padding
    if (i < s.size()) {
      res.push_back(to_fourbit_.at(s[i], NONE_CHAR));
    }

    return res;
  }

  std::string from_fourbit(const std::basic_string<char>& s, const size_t n)
  {
    // prepare the result string
    std::string res;
    res.reserve(n);

    // determine wether the packed string has padding
    const bool padded = (s.size() < n);

    // unpack
    for (size_t i = 0; i < s.size(); ++i) {
      auto char_pair = unpack_(s[i]);
      res.push_back(NT_MAP[char_pair.first]);
      if ((i < s.size() - 1u) or not padded) {
        res.push_back(NT_MAP[char_pair.second]);
      }
    }

    return res;
  }
  
private:
  Matrix<char> to_fourbit_; 
};