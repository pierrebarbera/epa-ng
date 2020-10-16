#include "Epatest.hpp"

#include "core/pll/rtree_mapper.hpp"
#include "io/file_io.hpp"
#include "seq/MSA.hpp"

using namespace std;

TEST( file_io, build_MSA_from_file )
{
  auto msa = build_MSA_from_file( env->reference_file, MSA_Info( env->reference_file ), true );

  EXPECT_EQ( msa.size(), 8ul );
  EXPECT_EQ( msa.num_sites(), 705ul );
}

TEST( file_io, file_check )
{
  EXPECT_ANY_THROW( file_check( "asjbjibvi.hhs" ) );

  EXPECT_NO_THROW( file_check( env->combined_file ) );
  EXPECT_NO_THROW( file_check( env->query_file ) );
  EXPECT_NO_THROW( file_check( env->reference_file ) );
  EXPECT_NO_THROW( file_check( env->tree_file ) );
}
