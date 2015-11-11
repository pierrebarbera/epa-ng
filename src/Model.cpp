#include "Model.h"

#include <vector>
#include <stdexcept>

using namespace std;

Model::Model(vector<double> base_frequencies, vector<double> substitution_rates)
  : base_frequencies_(base_frequencies), substitution_rates_(substitution_rates)
{

  if (this->base_frequencies().size() != STATES)
    throw runtime_error{"Inappropriate number of base frequencies"};

  //TODO what is the appropriate number of substitution rates?

}

Model::~Model()
{
  // vectors are automatically deallocated
}
