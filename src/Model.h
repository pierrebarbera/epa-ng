#include <vector>

//TODO find a better place for these. possible candidates for compile flag variation: DNA vs protein
#define STATES    4
#define RATE_CATS 4

/* Encapsulates the evolutionary model parameters
  TODO possible basepoint of model class hierarchy */
class Model
{
public:
	Model(std::vector<double> base_frequencies, std::vector<double> substitution_rates);
	~Model();

  std::vector<double> base_frequencies;
  std::vector<double> substitution_rates;

};