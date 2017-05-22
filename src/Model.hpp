#pragma once

#include <vector>
#include <unordered_map>
#include "constants.hpp"
#include "pll.h"

/* Encapsulates the evolutionary model parameters
  TODO possible basepoint of model class hierarchy */
class Model {
public:
  Model() : Model("DNA", "GTR", "EMPIRICAL") { }
  Model(std::string sequence_type, std::string model_id, std::string sub_matrix);
  ~Model() = default;

  // getters
  const std::vector<double>& base_frequencies() const {return base_frequencies_;}
  std::vector<int>& symmetries()  {return subs_symmetries_;}
  const std::vector<double>& substitution_rates() const {return substitution_rates_;}
  double alpha() const {return alpha_;}
  unsigned int states() const {return states_;}
  unsigned int rate_cats() const {return rate_cats_;}
  unsigned int const * char_map() const {return char_map_;}

  // setters
  void set_from_partition(pll_partition_t const * const partition);
  void base_frequencies(double* source, unsigned int length);
  void base_frequencies(std::vector<double> freqs);
  void substitution_rates(double* source, unsigned int length);
  void substitution_rates(std::vector<double> rates);
  void symmetries(int* source, unsigned int length);
  void alpha(double a) {alpha_ = a;};

private:
  unsigned int const * char_map_;
  unsigned int states_;
  unsigned int rate_cats_;
  double alpha_;
  std::vector<double> base_frequencies_;
  std::vector<double> substitution_rates_;
  /* six integer code defining the symmetries in the subsitution matrix
    examples:
    JC/F81:   000000
    HKY/K80:  010010
    */
  std::vector<int> subs_symmetries_;
};

unsigned int const * get_char_map(pll_partition_t const * const partition);
