#include "Epatest.hpp"

#include "src/set_manipulators.hpp"
#include "src/jplace_util.hpp"
#include "src/Range.hpp"

#include <vector>
#include <iostream>

using namespace std;

TEST(set_manipulators, split_sample)
{
  Sample sample_1;
  unsigned int s_a = 0, s_b = 1, s_c = 2;
  sample_1.emplace_back(s_a, 0);
  sample_1.back().emplace_back(1,-10,0.9,0.9);
  sample_1.emplace_back(s_b, 0);
  sample_1.back().emplace_back(2,-10,0.9,0.9);
  sample_1.emplace_back(s_c, 0);
  sample_1.back().emplace_back(3,-10,0.9,0.9);

  vector<Sample> parts;
  vector<vector<unsigned int>> split_map;
  split_map.push_back({s_a, s_b});
  split_map.push_back({s_c});

  assert(split_map[0].size() == 2);
  assert(split_map[1].size() == 1);

  split(sample_1, parts, split_map);

  ASSERT_EQ(2, parts.size());
  ASSERT_EQ(0, sample_1.size());

  ASSERT_EQ(2, parts[0].size());
  ASSERT_EQ(1, parts[1].size());

  EXPECT_EQ(1, parts[0][0][0].branch_id());
  EXPECT_EQ(2, parts[0][1][0].branch_id());
  EXPECT_EQ(3, parts[1][0][0].branch_id());
}

TEST(set_manipulators, split_sample_equal)
{
  Sample sample_1;
  unsigned int s_a = 0, s_b = 1, s_c = 2;
  sample_1.emplace_back(s_a, 0);
  sample_1.back().emplace_back(1,-10,0.9,0.9);
  sample_1.emplace_back(s_b, 0);
  sample_1.back().emplace_back(2,-10,0.9,0.9);
  sample_1.emplace_back(s_c, 0);
  sample_1.back().emplace_back(3,-10,0.9,0.9);

  vector<Sample> parts;

  split(sample_1, parts, 2);

  ASSERT_EQ(2, parts.size());
  ASSERT_NE(0, sample_1.size());

  ASSERT_EQ(2, parts[0].size());
  ASSERT_EQ(1, parts[1].size());

  EXPECT_EQ(1, parts[0][0][0].branch_id());
  EXPECT_EQ(2, parts[0][1][0].branch_id());
  EXPECT_EQ(3, parts[1][0][0].branch_id());
}

TEST(set_manipulators, split_work_equal)
{
  Sample sample_1;
  unsigned int s_a = 0, s_b = 1, s_c = 2;
  sample_1.emplace_back(s_a);
  sample_1.back().emplace_back(1,-10,0.9,0.9);
  sample_1.back().emplace_back(2,-10,0.9,0.9);
  sample_1.back().emplace_back(3,-10,0.9,0.9);
  sample_1.emplace_back(s_b);
  sample_1.back().emplace_back(2,-10,0.9,0.9);
  sample_1.emplace_back(s_c);
  sample_1.back().emplace_back(1,-10,0.9,0.9);
  sample_1.back().emplace_back(2,-10,0.9,0.9);
  sample_1.back().emplace_back(3,-10,0.9,0.9);

  // for (auto& pq : sample_1)
  // {
  //   printf("Sequence %d: ", pq.sequence_id());
  //   for (auto& p : pq)
  //     printf(" %d ", p.branch_id());
  //   printf("\n");
  // }

  // printf("\nWork");

  vector<Work> parts;
  Work work(sample_1);

  // for (auto i = work.begin(); i != work.end(); ++i)
  // {
  //   printf("\nbranch %d: ", i->first);
  //   for (auto& seq_id : i->second)
  //   {
  //     printf(" %d ", seq_id);
  //   }
  // }
  // printf("\n");

  split(work, parts, 2);
  
  auto total_placements = 0;
  for (auto& pq : sample_1)
    total_placements+=pq.size();

  auto total_work = 0;
  for (auto& w : parts)
    total_work+=w.size();

  EXPECT_EQ(total_placements, total_work);

  ASSERT_EQ(2, parts.size());
  ASSERT_EQ(total_placements, work.size());

  // for (auto& p : parts)
  // {
  //   printf("part:");
  //   for (auto i = p.begin(); i != p.end(); ++i)
  //   {
  //     printf("\nbranch %d: ", i->first);
  //     for (auto& seq_id : i->second)
  //     {
  //       printf(" %d ", seq_id);
  //     }
  //   }
  //   printf("\n");
  // }

  ASSERT_EQ(4, parts[0].size());
  ASSERT_EQ(3, parts[1].size());
}

TEST(set_manipulators, merge_work)
{
  Sample sample;
  unsigned int s_a = 0, s_b = 1, s_c = 2;
  sample.emplace_back(s_a);
  sample.back().emplace_back(1,-10,0.9,0.9);
  sample.back().emplace_back(2,-10,0.9,0.9);
  sample.back().emplace_back(3,-10,0.9,0.9);
  sample.emplace_back(s_b);
  sample.back().emplace_back(2,-10,0.9,0.9);
  sample.emplace_back(s_c);
  sample.back().emplace_back(1,-10,0.9,0.9);
  sample.back().emplace_back(2,-10,0.9,0.9);
  sample.back().emplace_back(3,-10,0.9,0.9);

  // printf("\nSample");

  // for (auto& pq : sample)
  // {
  //   printf("Sequence %d: ", pq.sequence_id());
  //   for (auto& p : pq)
  //     printf(" %d ", p.branch_id());
  //   printf("\n");
  // }

  vector<Work> parts;
  Work work(sample);
  split(work, parts, 2);

  // printf("\nWork");
  // for (auto i = work.begin(); i != work.end(); ++i)
  // {
  //   printf("\nbranch %d: ", i->first);
  //   for (auto& seq_id : i->second)
  //   {
  //     printf(" %d ", seq_id);
  //   }
  // }
  // printf("\n");


  Work dest;

  for(auto& w : parts) 
    merge(dest, w);

  // printf("\nDest");
  // for (auto i = dest.begin(); i != dest.end(); ++i)
  // {
  //   printf("\nbranch %d: ", i->first);
  //   for (auto& seq_id : i->second)
  //   {
  //     printf(" %d ", seq_id);
  //   }
  // }
  // printf("\n");

  EXPECT_EQ(work.size(), dest.size());

  for(auto& i : work) 
  {
    EXPECT_EQ(i.second.size(), dest[i.first].size());
  }


}

TEST(set_manipulators, merge_sample)
{
  // setup
  Sample sample_1;
  Sample sample_2;
  unsigned int s_a = 0, s_b = 1, s_c = 2, s_d = 3;
  sample_1.emplace_back(s_a, 0);
  sample_1.back().emplace_back(1,-10,0.9,0.9);
  sample_1.emplace_back(s_b, 0);
  sample_1.back().emplace_back(2,-10,0.9,0.9);
  sample_1.emplace_back(s_c, 0);
  sample_1.back().emplace_back(3,-10,0.9,0.9);

  assert(sample_1.size() == 3);
  assert(sample_1[0].size() == 1);
  assert(sample_1[1].size() == 1);
  assert(sample_1[2].size() == 1);

  sample_2.emplace_back(s_c, 0);
  sample_2.back().emplace_back(1,-10,0.9,0.9);
  sample_2.emplace_back(s_b, 0);
  sample_2.back().emplace_back(0,-10,0.9,0.9);
  sample_2.emplace_back(s_d, 0);
  sample_2.back().emplace_back(3,-10,0.9,0.9);

  assert(sample_2.size() == 3);
  assert(sample_2[0].size() == 1);
  assert(sample_2[1].size() == 1);
  assert(sample_2[2].size() == 1);

  merge(sample_1, sample_2);

  ASSERT_EQ(4, sample_1.size());

  EXPECT_EQ(1, sample_1[0].size());
  EXPECT_EQ(2, sample_1[1].size());
  EXPECT_EQ(2, sample_1[2].size());
  EXPECT_EQ(1, sample_1[3].size());

}

TEST(set_manipulators, find_collapse_equal_sequences)
{
  // buildup
  MSA msa;
  msa.append(string("T"),   string("AGCTAGCT"));
  msa.append(string("C"),   string("AGCCAGCT"));
  msa.append(string("G"),   string("AGCGAGCT"));
  msa.append(string("A"),   string("AGCAAGCT"));
  msa.append(string("t1"),  string("AGCTAGCT"));
  msa.append(string("c1"),  string("AGCCAGCT"));
  msa.append(string("g1"),  string("AGCGAGCT"));
  msa.append(string("a1"),  string("AGCAAGCT"));
  msa.append(string("t2"),  string("AGCTAGCT"));
  msa.append(string("c2"),  string("AGCCAGCT"));

  // test
  find_collapse_equal_sequences(msa);
  EXPECT_EQ(4, msa.size());

  EXPECT_EQ(3, msa[0].header_list().size());
  EXPECT_EQ(3, msa[1].header_list().size());
  EXPECT_EQ(2, msa[2].header_list().size());
  EXPECT_EQ(2, msa[3].header_list().size());
}

TEST(set_manipulators, get_valid_range)
{
  string s1("---------GGGCCCGTAT-------");//(9,19)
  string s2("GGGCCCGTAT-------");         //(0,10)
  string s3("-GGGC---CCG-TAT");           //(1,15)

  Range r;
  r = get_valid_range(s1);
  EXPECT_EQ(9, r.begin);
  EXPECT_EQ(10, r.span);

  r = get_valid_range(s2);
  EXPECT_EQ(0, r.begin);
  EXPECT_EQ(10, r.span);

  r = get_valid_range(s3);
  EXPECT_EQ(1, r.begin);
  EXPECT_EQ(14, r.span);
}

TEST(set_manipulators, discard_bottom_x_percent)
{
  // setup
  Sample sample;
  unsigned int s_a = 0, s_b = 1, s_c = 2;
  sample.emplace_back(s_a, 0);
  vector<double> weights_a({0.001,0.23,0.05,0.02,0.4,0.009,0.2,0.09});
  vector<double> weights_b({0.01,0.02,0.005,0.002,0.94,0.003,0.02});
  unsigned int num_expected[3] = {4,4,1};

  for (auto n : weights_a) {
    sample.back().emplace_back(1,-10,0.9,0.9);
    sample.back().back().lwr(n);
  }
  sample.emplace_back(s_b, 0);
  for (auto n : weights_b) {
    sample.back().emplace_back(1,-10,0.9,0.9);
    sample.back().back().lwr(n);
  }
  sample.emplace_back(s_c, 0);
  sample.back().emplace_back(1,-10,0.9,0.9);
  sample.back().back().lwr(1.0);

  // tests
  discard_bottom_x_percent(sample, 0.5);
  int i =0;
  for (auto pq : sample) {
    unsigned int num = 0;
    for (auto p : pq) {
      (void)p;
      num++;
    }
    EXPECT_EQ(num_expected[i++], num);
  }
}

TEST(set_manipulators, discard_by_accumulated_threshold)
{
  // setup
  Sample sample;
  unsigned int s_a = 0, s_b = 1, s_c = 2;
  sample.emplace_back(s_a, 0);
  vector<double> weights_a({0.001,0.23,0.05,0.02,0.4,0.009,0.2,0.09});
  vector<double> weights_b({0.01,0.02,0.005,0.002,0.94,0.003,0.02});
  unsigned int num_expected[3] = {5,2,1};

  for (auto n : weights_a) {
    sample.back().emplace_back(1,-10,0.9,0.9);
    sample.back().back().lwr(n);
  }
  sample.emplace_back(s_b, 0);
  for (auto n : weights_b) {
    sample.back().emplace_back(1,-10,0.9,0.9);
    sample.back().back().lwr(n);
  }
  sample.emplace_back(s_c, 0);
  sample.back().emplace_back(1,-10,0.9,0.9);
  sample.back().back().lwr(1.0);

  // tests
  discard_by_accumulated_threshold(sample, 0.95);
  int i =0;
  for (auto pq : sample) {
    unsigned int num = 0;
    for (auto p : pq) {
      (void)p;
      num++;
    }
    EXPECT_EQ(num_expected[i++], num);
  }
}

TEST(set_manipulators, discard_by_support_threshold)
{
  // setup
  Sample sample;
  unsigned int s_a = 0, s_b = 1, s_c = 2;
  sample.emplace_back(s_a, 0);
  vector<double> weights_a{0.001,0.23,0.05,0.02,0.4,0.009,0.2,0.09};
  vector<double> weights_b{0.01,0.02,0.005,0.002,0.94,0.003,0.02};
  unsigned int num_expected[3] = {6,3,1};

  for (auto n : weights_a) {
    sample.back().emplace_back(1,-10,0.9,0.9);
    sample.back().back().lwr(n);
  }
  sample.emplace_back(s_b, 0);
  for (auto n : weights_b) {
    sample.back().emplace_back(1,-10,0.9,0.9);
    sample.back().back().lwr(n);
  }
  sample.emplace_back(s_c, 0);
  sample.back().emplace_back(1,-10,0.9,0.9);
  sample.back().back().lwr(1.0);

  // tests
  discard_by_support_threshold(sample, 0.01);
  int i =0;
  for (auto pq : sample) {
    unsigned int num = 0;
    for (auto p : pq) {
      (void)p;
      num++;
    }
    EXPECT_EQ( num_expected[i++], num);
  }
}
