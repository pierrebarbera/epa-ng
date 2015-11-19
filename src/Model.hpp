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
	Model(std::vector<double> base_frequencies, std::vector<double> substitution_rates, double alpha);
	~Model();
  inline const std::vector<double>& base_frequencies() const {return base_frequencies_;};
  inline const std::vector<double>& substitution_rates() const {return substitution_rates_;};
  inline double alpha() const {return alpha_;};

private:
  double alpha_;
  std::vector<double> base_frequencies_;
  std::vector<double> substitution_rates_;

};

#endif
