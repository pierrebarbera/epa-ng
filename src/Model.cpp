#include "Model.hpp"

#include <vector>
#include <stdexcept>

using namespace std;

Model::Model(vector<double> base_frequencies, vector<double> substitution_rates,
  double alpha, std::vector<int> symmetries)
  : alpha_(alpha), base_frequencies_(base_frequencies), substitution_rates_(substitution_rates),
    subs_symmetries_(symmetries)
{

  if (base_frequencies_.size() != STATES)
    throw runtime_error{"Inappropriate number of base frequencies"};

  // TODO what is the appropriate number of substitution rates?
  // TODO what is the apropriate alpha input range

}

void Model::base_frequencies(double *source, unsigned int length)
{
  if (base_frequencies_.size() != length)
    throw runtime_error{"Inappropriate number of base frequencies"};

  for (unsigned int i = 0; i < length; ++i)
    base_frequencies_[i] = source[i];

}

void Model::substitution_rates(double *source, unsigned int length)
{
  if (substitution_rates_.size() != length)
    throw runtime_error{"Inappropriate number of substitution rates"};

  for (unsigned int i = 0; i < length; ++i)
    substitution_rates_[i] = source[i];

}
