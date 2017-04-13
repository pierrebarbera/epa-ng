#include "Model.hpp"

#include <stdexcept>
#include <string>
#include <algorithm>

#include "pllhead.hpp"

// map for determining model symmetries
static const std::unordered_map<std::string, std::vector<int>> MODEL_MAP(
  {
    {"JC69", {0,0,0,0,0,0}},
    {"K80", {0,1,0,0,1,0}},
    {"GTR", {0,1,2,3,4,5}}
  }
);

// map for determining rates and frequencies symmetries
static const std::unordered_map<std::string, std::pair<const double*, const double*>> AA_RATE_FREQ_MAP(
  {
    {"DAYHOFF", {pll_aa_rates_dayhoff, pll_aa_freqs_dayhoff}},
    {"LG", {pll_aa_rates_lg, pll_aa_freqs_lg}},
    {"DCMUT", {pll_aa_rates_dcmut, pll_aa_freqs_dcmut}},
    {"JTT", {pll_aa_rates_jtt, pll_aa_freqs_jtt}},
    {"MTREV", {pll_aa_rates_mtrev, pll_aa_freqs_mtrev}},
    {"WAG", {pll_aa_rates_wag, pll_aa_freqs_wag}},
    {"RTREV", {pll_aa_rates_rtrev, pll_aa_freqs_rtrev}},
    {"CPREV", {pll_aa_rates_cprev, pll_aa_freqs_cprev}},
    {"VT", {pll_aa_rates_vt, pll_aa_freqs_vt}},
    {"BLOSUM62", {pll_aa_rates_blosum62, pll_aa_freqs_blosum62}},
    {"MTMAM", {pll_aa_rates_mtmam, pll_aa_freqs_mtmam}},
    {"MTART", {pll_aa_rates_mtart, pll_aa_freqs_mtart}},
    {"MTZOA", {pll_aa_rates_mtzoa, pll_aa_freqs_mtzoa}},
    {"PMB", {pll_aa_rates_pmb, pll_aa_freqs_pmb}},
    {"HIVB", {pll_aa_rates_hivb, pll_aa_freqs_hivb}},
    {"HIVW", {pll_aa_rates_hivw, pll_aa_freqs_hivw}},
    {"JTTDCMUT", {pll_aa_rates_jttdcmut, pll_aa_freqs_jttdcmut}},
    {"FLU", {pll_aa_rates_flu, pll_aa_freqs_flu}},
    {"STMTREV", {pll_aa_rates_stmtrev, pll_aa_freqs_stmtrev}}
  }
);

Model::Model(std::string sequence_type, std::string model_id, std::string sub_matrix)
  : alpha_(1.0), base_frequencies_({0.25,0.25,0.25,0.25}),
    substitution_rates_({0.5,0.5,0.5,0.5,0.5,1.0})
{
  
  // tolerate case insensitivity
  transform(model_id.begin(), model_id.end(),model_id.begin(), ::toupper);
  transform(sequence_type.begin(), sequence_type.end(),sequence_type.begin(), ::toupper);
  transform(sub_matrix.begin(), sub_matrix.end(),sub_matrix.begin(), ::toupper);

  if(sequence_type.compare("DNA") != 0 
    and sequence_type.compare("AA") != 0)
    throw std::runtime_error{std::string("Sequence data type not recognized! Input: ") + sequence_type};


  if (sequence_type.compare("DNA") == 0)
  {
    states_= 4;
    rate_cats_ = 4;
    char_map_ = pll_map_nt;
  }
  else if (sequence_type.compare("AA") == 0)
  {
    states_ = 20;
    rate_cats_ = 4;
    char_map_ = pll_map_aa;

    auto find_iter = AA_RATE_FREQ_MAP.find(sub_matrix);
    if (find_iter == AA_RATE_FREQ_MAP.end())
      throw std::runtime_error{std::string("Sub. Matrix Identifier not found! String passed: ") + sub_matrix};

    auto rates = find_iter->second.first;
    auto freqs = find_iter->second.second;

    base_frequencies_.clear();
    for (size_t i = 0; i < 20; ++i)
      base_frequencies_.push_back(freqs[i]);
    for (size_t i = 0; i < 190; ++i)
      substitution_rates_.push_back(rates[i]);

  }
  
  // set symmetries
  auto find_iter = MODEL_MAP.find(model_id);
  if (find_iter == MODEL_MAP.end())
    throw std::runtime_error{std::string("Model Identifier not found! String passed: ") + model_id};
  for (auto n : find_iter->second)
    subs_symmetries_.push_back(n);
}

void Model::base_frequencies(double *source, unsigned int length)
{
  if (base_frequencies_.size() != length) {
    throw std::runtime_error{"Inappropriate number of base frequencies"};
  }

  for (size_t i = 0; i < length; ++i) {
    base_frequencies_[i] = source[i];
  }
}

void Model::base_frequencies(std::vector<double> freqs)
{
  base_frequencies(&freqs[0], freqs.size());
}

void Model::substitution_rates(double *source, unsigned int length)
{
  if (substitution_rates_.size() != length) {
    throw std::runtime_error{"Inappropriate number of substitution rates"};
  }

  for (size_t i = 0; i < length; ++i) {
    substitution_rates_[i] = source[i];
  }
}

void Model::substitution_rates(std::vector<double> rates)
{
  substitution_rates(&rates[0], rates.size());
}

void Model::symmetries(int* source, unsigned int length)
{
  if (subs_symmetries_.size() != length)
    throw std::runtime_error{"Inappropriate number of substitution symmetries"};

  for (size_t i = 0; i < length; ++i) {
    subs_symmetries_[i] = source[i];
  }
}
