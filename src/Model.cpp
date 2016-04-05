#include "Model.hpp"

#include <stdexcept>
#include <string>
#include <algorithm>

using namespace std;

Model::Model(string model_id)
  : alpha_(1.0), base_frequencies_({0.25,0.25,0.25,0.25}),
  substitution_rates_({0.5,0.5,0.5,0.5,0.5,1.0})
{
  // tolerate case insensitivity
  transform(model_id.begin(), model_id.end(),model_id.begin(), ::toupper);
  // set symmetries
  // subs_symmetries_ = MODEL_MAP[model_id.c_str()];
  auto find_iter = MODEL_MAP.find(model_id);
  if (find_iter == MODEL_MAP.end())
    throw runtime_error{string("Model Identifier not found! String passed: ") + model_id};
  for (auto n : find_iter->second)
    subs_symmetries_.push_back(n);
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

void Model::symmetries(int* source, unsigned int length)
{
  if (subs_symmetries_.size() != length)
    throw runtime_error{"Inappropriate number of substitution symmetries"};

  for (unsigned int i = 0; i < length; ++i)
    subs_symmetries_[i] = source[i];

}
