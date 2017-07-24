#include "Epatest.hpp"

#include "src/Binary_Fasta.hpp"
#include "src/file_io.hpp"
#include "src/MSA.hpp"

#include "genesis/utils/core/options.hpp"

static void compare_msas(const MSA& lhs, const MSA& rhs)
{
  ASSERT_EQ(lhs.size(), rhs.size());

  for (size_t i = 0; i < lhs.size(); ++i) {
    EXPECT_STREQ(lhs[i].header().c_str(), rhs[i].header().c_str());
    EXPECT_STREQ(lhs[i].sequence().c_str(), rhs[i].sequence().c_str());
  }
}

TEST(Binary_Fasta, 4bit_store_and_load)
{
  genesis::utils::Options::get().allow_file_overwriting(true);

  Binary_Fasta conv;

  auto msa = build_MSA_from_file(env->query_file);

  std::string binfile_name = env->query_file + ".bin";

  conv.save(msa, binfile_name);

  auto read_msa = conv.load(binfile_name);

  compare_msas(msa, read_msa);
}