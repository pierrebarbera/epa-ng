#include "Epatest.hpp"
#include "check_equal.hpp"
#include "print.hpp"
#include "get_randomized.hpp"

#include "core/BranchBuffer.hpp"

#include "core/pll/pllhead.hpp"
#include "core/Work.hpp"
#include "seq/MSA.hpp"
#include "seq/MSA_Info.hpp"
#include "tree/Tree.hpp"
#include "tree/Tiny_Tree.hpp"
#include "util/Options.hpp"
#include "util/parse_model.hpp"
#include "util/memory.hpp"
#include "io/file_io.hpp"


static void test_buffer_impl( std::string const& tree_file,
                              std::string const& msa_file,
                              std::string const& model_file,
                              Options const options,
                              bool const whitelisted )
{
  Work work;

  // buildup
  auto ref_info = make_msa_info( msa_file );
  MSA_Info qry_info = ref_info;
  qry_info.sequences( 0 );
  auto model      = raxml::Model( parse_model( model_file ) );
  auto footprint  = Memory_Footprint( ref_info, qry_info, model, options );
  auto msa        = build_MSA_from_file( msa_file,
                                  ref_info,
                                  options.premasking );
  auto ref_tree   = Tree( tree_file, msa, model, options, footprint );

  Options def_opts;
  auto true_tree  = Tree( tree_file, msa, model, def_opts, footprint );
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

  size_t const block_size = ref_tree.memsave().concurrent_branches;
  BranchBuffer bufferino( &ref_tree, block_size, work );

  ASSERT_EQ( std::accumulate( std::begin( bufferino.whitelist() ),
                              std::end( bufferino.whitelist() ),
                              0ul ),
             work.branches() );

  BranchBuffer::container_type block;

  std::set< unsigned int > branches_tested;
  size_t sum_branches = 0;

  while( bufferino.get_next( block ) ) {
    for( auto& buffered_tt : block ) {
      sum_branches++;
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

      branches_tested.insert( branch_id );
    }
  }

  // check we tested the correct number of branches
  ASSERT_EQ( sum_branches,
             whitelisted ? work.branches() : ref_tree.nums().branches );
  ASSERT_EQ( branches_tested.size(),
             whitelisted ? work.branches() : ref_tree.nums().branches );
}

static void test_buffer( Options const options,
                         std::string const& tree_file  = env->tree_file,
                         std::string const& msa_file   = env->reference_file,
                         std::string const& model_file = env->model_file )
{
  test_buffer_impl( tree_file, msa_file, model_file, options, false );
}

static void
test_buffer_whitelisted( Options const options,
                         std::string const& tree_file  = env->tree_file,
                         std::string const& msa_file   = env->reference_file,
                         std::string const& model_file = env->model_file )

{
  test_buffer_impl( tree_file, msa_file, model_file, options, true );
}

TEST( BranchBuffer, buffer )
{
  Options options;

  options.memsave = Memsave_Option::Mode::kFull;

  test_buffer( options );
}

TEST( BranchBuffer, buffer_whitelisted )
{
  Options options;

  options.memsave = Memsave_Option::Mode::kFull;

  test_buffer_whitelisted( options );
}

TEST( BranchBuffer, neotrop )
{
  Options options;
  options.memsave = Memsave_Option::Mode::kFull;

  std::string dir = env->data_dir + "neotrop/";
  test_buffer( options,
               dir + "tree.newick",
               dir + "reference.fasta.gz",
               dir + "infofile" );
}

TEST( BranchBuffer, neotrop_whitelisted )
{
  Options options;
  options.memsave = Memsave_Option::Mode::kFull;

  std::string dir = env->data_dir + "neotrop/";
  test_buffer_whitelisted( options,
                           dir + "tree.newick",
                           dir + "reference.fasta.gz",
                           dir + "infofile" );
}
