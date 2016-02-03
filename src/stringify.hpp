#pragma once

#include <string>
#include <sstream>

#include "Model.hpp"

constexpr char NEWL = '\n';

void rwnd(std::ostringstream& ss, unsigned int chars);

std::string to_string(Model& model);
