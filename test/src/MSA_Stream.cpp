#include "Epatest.hpp"

#include "seq/MSA_Stream.hpp"
#include "seq/MSA.hpp"
#include "io/file_io.hpp"

#include <string>

using namespace std;

TEST(MSA_Stream, reading)
{
  MSA complete_msa = build_MSA_from_file(env->combined_file);
  const auto chunk_size = 3;
  MSA read_msa;
  MSA_Stream streamed_msa(env->combined_file, chunk_size);

  for (size_t i = 0; i < complete_msa.size(); i++)
  {
    if ((i % chunk_size) == 0)
    {
      streamed_msa.read_next(read_msa, chunk_size);
    }
    EXPECT_EQ(complete_msa[i], read_msa[i % chunk_size]);
  }
  MSA_Stream dummy;
}
