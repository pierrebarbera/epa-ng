#ifndef EPA_JPLACE_UTIL_H_
#define EPA_JPLACE_UTIL_H_

#include <string>
#include <sstream>

#include "PQuery.hpp"
#include "PQuery_Set.hpp"

#define TAB "  ";
#define NEWL "\n";

std::string pquery_to_jplace_string(const PQuery& p);
std::string pquery_set_to_jplace_string(const PQuery_Set& ps, std::string& invocation);

#endif
