#include "Model.hpp"

#include <vector>
#include <stdexcept>

using namespace std;

Model::Model(vector<double> base_frequencies, vector<double> substitution_rates,
  double alpha, std::vector<int> symmetries)
  : alpha_(alpha), base_frequencies_(base_frequencies), substitution_rates_(substitution_rates),
    subs_symmetries_(symmetries)
{

  if (this->base_frequencies().size() != STATES)
    throw runtime_error{"Inappropriate number of base frequencies"};

  // TODO what is the appropriate number of substitution rates?
  // TODO what is the apropriate alpha input range

}

Model::~Model()
{
  // vectors are automatically deallocated
}
