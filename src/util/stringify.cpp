#include "util/stringify.hpp"

std::string stringify(raxml::Model& model)
{
  std::ostringstream output;

  output << "Substitution Matrix Symmetries: " << NEWL;
  output << stringify(model.submodel(0).rate_sym()) << NEWL;

  output << "Base Frequencies: " << NEWL;
  output << stringify(model.base_freqs(0)) << NEWL;

  output << "Substitution Rates: " << NEWL;
  output << stringify(model.subst_rates(0)) << NEWL;

  output << "Alpha: " << model.alpha() << NEWL;

  return output.str();
}
