#pragma once

#include <stdexcept>
#include <string>

#include "core/pll/pllhead.hpp"

inline void handle_pll_failure( bool const failed, std::string const message )
{
  if( failed ) {
    throw std::runtime_error{ message + " (pll says: " + pll_errmsg + ")" };
  }
}