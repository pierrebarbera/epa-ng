#pragma once

#include <string>
#include <fstream>
#include <sstream>
#include <vector>

#include "Model.hpp"

constexpr char NEWL = '\n';

template <class T, class lambda>
std::string stringify_vec_impl( const std::vector<T>& vec, 
                                lambda toString)
{
  std::ostringstream output;

  size_t i = 0;
  for (auto& e : vec) {
    output << toString(e);
    if (++i < vec.size()) {
      output << ", ";  
    }
  }

  return output.str();
}

std::string stringify(Model& model);

template <class T>
std::string stringify(const std::vector<T>& vec)
{
  return stringify_vec_impl(vec, [](const auto& elem) {
    return elem;
  });
}

template <class T>
std::string stringify(const std::vector<std::vector<T>>& vec)
{
  return stringify_vec_impl(vec, [](const auto& elem) {
    return elem.size();
  });
}
