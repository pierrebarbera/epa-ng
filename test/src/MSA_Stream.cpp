#include "Epatest.hpp"

#include "src/MSA_Stream.hpp"
#include "src/MSA.hpp"
#include "src/file_io.hpp"

#include <string>

using namespace std;

TEST(MSA_Stream, reading)
{
  MSA complete_msa = build_MSA_from_file(env->combined_file);
  MSA_Stream streamed_msa(env->combined_file);
  auto chunk_size = 3;

  for (size_t i = 0; i < complete_msa.size(); i++)
  {
    if ((i % chunk_size) == 0)
    {
      streamed_msa.clear();
      streamed_msa.read_next(chunk_size);
    }
    EXPECT_EQ(complete_msa[i], streamed_msa[i % chunk_size]);
  }
}
