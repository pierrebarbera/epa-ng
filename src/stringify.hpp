#pragma once

#include <string>
#include <fstream>
#include <sstream>

#include "Model.hpp"

constexpr char NEWL = '\n';

void rwnd(std::ostream& ss, unsigned int chars);

std::string to_string(Model& model);
