#include "Epatest.hpp"

#include "src/calculation.hpp"
#include "src/jplace_util.hpp"

#include <vector>
#include <iostream>

using namespace std;

TEST(calculation, discard_by_accumulated_threshold)
{
  // setup
  PQuery_Set pqs;
  Sequence s_a, s_b, s_c;
  pqs.emplace_back(s_a, 0);
  vector<double> weights_a({0.001,0.23,0.05,0.02,0.4,0.009,0.2,0.09});
  vector<double> weights_b({0.01,0.02,0.005,0.002,0.94,0.003,0.02});
  unsigned int num_expected[3] = {5,2,1};

  for (auto n : weights_a) {
    pqs.back().emplace_back(1,-10,0.9,0.9);
    pqs.back().back().lwr(n);
  }
  pqs.emplace_back(s_b, 0);
  for (auto n : weights_b) {
    pqs.back().emplace_back(1,-10,0.9,0.9);
    pqs.back().back().lwr(n);
  }
  pqs.emplace_back(s_c, 0);
  pqs.back().emplace_back(1,-10,0.9,0.9);
  pqs.back().back().lwr(1.0);

  // tests
  discard_by_accumulated_threshold(pqs, 0.95);
  int i =0;
  for (auto pq : pqs) {
    unsigned int num = 0;
    for (auto p : pq) {
      (void)p;
      num++;
    }
    EXPECT_EQ(num, num_expected[i++]);
  }

  // string inv("blorp");
  // cout << pquery_set_to_jplace_string(pqs, inv);

}

TEST(calculation, discard_by_support_threshold)
{
  // setup
  PQuery_Set pqs;
  Sequence s_a, s_b, s_c;
  pqs.emplace_back(s_a, 0);
  vector<double> weights_a{0.001,0.23,0.05,0.02,0.4,0.009,0.2,0.09};
  vector<double> weights_b{0.01,0.02,0.005,0.002,0.94,0.003,0.02};
  unsigned int num_expected[3] = {6,3,1};

  for (auto n : weights_a) {
    pqs.back().emplace_back(1,-10,0.9,0.9);
    pqs.back().back().lwr(n);
  }
  pqs.emplace_back(s_b, 0);
  for (auto n : weights_b) {
    pqs.back().emplace_back(1,-10,0.9,0.9);
    pqs.back().back().lwr(n);
  }
  pqs.emplace_back(s_c, 0);
  pqs.back().emplace_back(1,-10,0.9,0.9);
  pqs.back().back().lwr(1.0);

  // tests
  discard_by_support_threshold(pqs, 0.01);
  int i =0;
  for (auto pq : pqs) {
    unsigned int num = 0;
    for (auto p : pq) {
      (void)p;
      num++;
    }
    EXPECT_EQ(num, num_expected[i++]);
  }

}
