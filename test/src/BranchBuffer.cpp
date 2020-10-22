#include "Epatest.hpp"
#include "check_equal.hpp"
#include "print.hpp"
#include "get_randomized.hpp"

#include "core/BranchBuffer.hpp"

#include "core/Work.hpp"
#include "seq/MSA.hpp"
#include "seq/MSA_Info.hpp"
#include "tree/Tree.hpp"
#include "tree/Tiny_Tree.hpp"
#include "util/Options.hpp"
#include "io/file_io.hpp"
#include "core/pll/pllhead.hpp"

static void test_buffer_impl( Options const options, bool const whitelisted )
{
  Work work;

  // buildup
  auto msa        = build_MSA_from_file( env->reference_file,
                                  MSA_Info( env->reference_file ),
                                  options.premasking );
  auto ref_tree   = Tree( env->tree_file, msa, env->model, options );

  Options def_opts;
  auto true_tree  = Tree( env->tree_file, msa, env->model, def_opts );
  std::vector< pll_unode_t* > branches( true_tree.nums().branches );
  auto num_traversed_branches
      = utree_query_branches( true_tree.tree(), &branches[ 0 ] );

  ASSERT_EQ( num_traversed_branches, true_tree.nums().branches );

  check_equal( *ref_tree.tree(), *true_tree.tree() );

  ASSERT_TRUE( ref_tree.branch_id() == true_tree.branch_id() );

  // create the buffer  
  if( whitelisted ) {
    work = get_randomized_Work( ref_tree.nums().branches, 1, 0.5 );
  }

  size_t const block_size = 10;
  BranchBuffer bufferino( &ref_tree, block_size, work );

  BranchBuffer::container_type block;

  size_t branches_tested = 0;

  while( bufferino.get_next( block ) ) {
    for( auto& buffered_tt : block ) {
      auto const branch_id = buffered_tt.branch_id();
      // check for signs of demonic corruption
      ASSERT_NE( branch_id, -1u );

      auto const true_tree_branch_id
          = true_tree.branch_id( branches[ branch_id ]->node_index );
      ASSERT_EQ( true_tree_branch_id, branch_id );

      // make a Tiny_Tree from the old-style ref_tree, for the branch_id of the
      // current Tiny Tree coming from the buffer
      Tiny_Tree true_tt(  branches[ branch_id ],
                          branch_id,
                          true_tree );

      // print_tipchars( *buffered_tt.partition() );
      // print_tipchars( *true_tt.partition() );

      // compare the two
      check_equal( buffered_tt, true_tt );

      branches_tested++;
    }
  }

  EXPECT_EQ( branches_tested,
             whitelisted ? work.branches() : ref_tree.nums().branches );
}

static void test_buffer( Options const options )
{
  test_buffer_impl( options, false );
}

static void test_buffer_whitelisted( Options const options )
{
  test_buffer_impl( options, true );
}

TEST( BranchBuffer, buffer )
{
  Options options;

  options.memsave = true;

  test_buffer( options );
  // all_combinations( test_buffer );
}

TEST( BranchBuffer, buffer_whitelisted )
{
  Options options;

  options.memsave = true;

  test_buffer_whitelisted( options );
  // all_combinations( test_buffer );
}
