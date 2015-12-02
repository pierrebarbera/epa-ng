#ifndef EPA_MODEL_H_
#define EPA_MODEL_H_

#include <vector>

//TODO find a better place for these. possible candidates for compile flag variation: DNA vs protein
#define STATES    4
#define RATE_CATS 4

/* Encapsulates the evolutionary model parameters
  TODO possible basepoint of model class hierarchy */
class Model
{
public:
  Model() = default;
	Model(std::vector<double> base_frequencies, std::vector<double> substitution_rates,
    double alpha, std::vector<int> symmetries = {0,0,0,0,0,0});
	~Model();
  inline const std::vector<double>& base_frequencies() const {return base_frequencies_;};
  inline std::vector<int>& symmetries() {return subs_symmetries_;};
  inline const std::vector<double>& substitution_rates() const {return substitution_rates_;};
  inline double alpha() const {return alpha_;};

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

#endif
