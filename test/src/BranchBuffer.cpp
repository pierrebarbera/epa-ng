#include "Epatest.hpp"

#include "core/BranchBuffer.hpp"

#include "seq/MSA.hpp"
#include "tree/Tree.hpp"

static void test_buffer( Options const options )
{
  // buildup
  auto msa      = build_MSA_from_file( env->reference_file,
                                  MSA_Info( env->reference_file ),
                                  options.premasking );
  auto ref_tree = Tree( env->tree_file, msa, env->model, options );

  // create the buffer
  size_t const block_size = 10;
  BranchBuffer bufferino( &ref_tree, block_size );

  BranchBuffer::container_type block;

  while( bufferino.get_next( block, block_size ) ) {
    std::cout << block.size() << std::endl;
  }
}

TEST( BranchBuffer, buffer )
{
  Options options;

  options.memsave = true;

  test_buffer( options );
  // all_combinations( test_buffer );
}
