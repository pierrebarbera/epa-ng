#include "stringify.hpp"

using namespace std;

std::string to_string(Model& model)
{
  ostringstream output;
  size_t i = 0;

  output << "Substitution Matrix Symmetries: " << NEWL;
  for (auto sym : model.symmetries()) {
    output << to_string(sym);
    if (++i < model.symmetries().size()) {
      output << ",";  
    }
  }

  output << NEWL;

  output << "Base Frequencies: " << NEWL;
  i = 0;
  for (auto freq : model.base_frequencies()) {
    output << to_string(freq);
    if (++i < model.base_frequencies().size()) {
      output << ",";  
    }
  }

  output << NEWL;

  output << "Substitution Rates: " << NEWL;
  i = 0;
  for (auto rate : model.substitution_rates()){
    output << to_string(rate) << ", ";
    if (++i < model.substitution_rates().size()) {
      output << ",";  
    }
  }

  output << NEWL;

  output << "Alpha: " << model.alpha() << NEWL;

  return output.str();
}
