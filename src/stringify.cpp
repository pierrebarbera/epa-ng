#include "stringify.hpp"

using namespace std;

void rwnd(ostringstream& ss, unsigned int chars)
{
  long pos = ss.tellp();
  ss.seekp (pos - (chars * sizeof(char)));
}

std::string to_string(Model& model)
{
  ostringstream output;

  output << "Substitution Matrix Symmetries: " << NEWL;
  for (auto sym : model.symmetries())
    output << to_string(sym) << ", ";

  rwnd(output, 2);
  output << NEWL;

  output << "Base Frequencies: " << NEWL;
  for (auto freq : model.base_frequencies())
    output << to_string(freq) << ", ";

  rwnd(output, 2);
  output << NEWL;

  output << "Substitution Rates: " << NEWL;
  for (auto rate : model.substitution_rates())
    output << to_string(rate) << ", ";

  rwnd(output, 2);
  output << NEWL;

  output << "Alpha: " << model.alpha() << NEWL;

  return output.str();
}
