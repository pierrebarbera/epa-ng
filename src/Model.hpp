#pragma once

#include <vector>
#include <unordered_map>

// TODO find a better place for these. possible candidates for compile flag variation: DNA vs protein
#define STATES    4
#define RATE_CATS 4

// map for determining model symmetries
const std::unordered_map<std::string, std::vector<int>> MODEL_MAP(
  {
    {"JC69", {0,0,0,0,0,0}},
    {"K80", {0,1,0,0,1,0}},
    {"GTR", {0,1,2,3,4,5}}
  }
);

/* Encapsulates the evolutionary model parameters
  TODO possible basepoint of model class hierarchy */
class Model
{
public:
  Model() = default;
  Model(std::string model_id);
	// Model(std::vector<double> base_frequencies, std::vector<double> substitution_rates,
  //   double alpha, std::vector<int> symmetries = {0,0,0,0,0,0});
	~Model() = default;

  // getters
  const std::vector<double>& base_frequencies() const {return base_frequencies_;};
  std::vector<int>& symmetries()  {return subs_symmetries_;};
  const std::vector<double>& substitution_rates() const {return substitution_rates_;};
  double alpha() const {return alpha_;};

  // setters
  void base_frequencies(double* source, unsigned int length);
  void substitution_rates(double* source, unsigned int length);
  void alpha(double a) {alpha_ = a;};

private:
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
