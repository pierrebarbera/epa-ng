#include "stdinc.h"

using namespace std;

//TODO find a better place for these. possible candidates for compile flag variation: DNA vs protein
#define STATES    4
#define RATE_CATS 4

/* Encapsulates the evolutionary model parameters
  TODO possible basepoint of model class hierarchy */
class Model
{
public:
	Model(vector<double> base_frequencies, vector<double> substitution_rates);
	~Model();

  vector<double> base_frequencies;
  vector<double> substitution_rates;

};