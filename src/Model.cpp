#include "Model.hpp"

#include <stdexcept>
#include <string>
#include <algorithm>

using namespace std;

// Model::Model(vector<double> base_frequencies, vector<double> substitution_rates,
//   double alpha, std::vector<int> symmetries)
//   : alpha_(alpha), base_frequencies_(base_frequencies), substitution_rates_(substitution_rates),
//     subs_symmetries_(symmetries)
// {
//
//   if (base_frequencies_.size() != STATES)
//     throw runtime_error{"Inappropriate number of base frequencies"};
//
//   // TODO what is the appropriate number of substitution rates?
//   // TODO what is the apropriate alpha input range
//
// }

Model::Model(string model_id)
  : alpha_(1.0), base_frequencies_({0.25,0.25,0.25,0.25}),
  substitution_rates_({0.5,0.5,0.5,0.5,0.5,1.0})
{
  // tolerate case insensitivity
  transform(model_id.begin(), model_id.end(),model_id.begin(), ::toupper);
  // set symmetries
  auto find_iter = model_map_.find(model_id);
  if (find_iter == model_map_.end())
    throw runtime_error{string("Model Identifier not found! String passed: ") + model_id};
  for (auto n : find_iter->second)
    subs_symmetries_.push_back(n);
}

Model::model_map_t Model::model_map_ =
{
  {"JC69", {0,0,0,0,0,0}},
  {"K80", {0,1,0,0,1,0}},
  {"GTR", {0,1,2,3,4,5}}
};

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
