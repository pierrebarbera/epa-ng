#include "Epatest.hpp"

#include "check_equal.hpp"

#include "core/place.hpp"
#include "core/raxml/Model.hpp"
#include "io/Binary_Fasta.hpp"
#include "io/file_io.hpp"
#include "seq/MSA.hpp"
#include "seq/MSA_Info.hpp"
#include "tree/Tree.hpp"
#include "util/Options.hpp"
#include "util/parse_model.hpp"

#include "genesis/genesis.hpp"

#include <limits>
#include <string>
#include <vector>

using namespace genesis;

/**
 * Set of tests performing regression testing: for each tested input, a
 * previously calculated output exists. We re-run the specified test, then check
 * if the results have changed.
 */

static void full_test( std::string const& previous_jplace,
                       std::string const& tree_file,
                       std::string const& msa_file,
                       std::string const& query_file,
                       std::string const& model_file,
                       Options const& options )
{
  utils::Options::get().allow_file_overwriting( true );

  std::string invocation( "./this --is -a test" );
  std::string outdir = env->out_dir;
  MSA_Info qry_info( query_file );
  MSA_Info ref_info( msa_file );
  MSA_Info::or_mask( ref_info, qry_info );

  auto msa = build_MSA_from_file( msa_file, ref_info, options.premasking );

  auto model = raxml::Model( parse_model( model_file ) );

  Tree ref_tree( tree_file, msa, model, options );

  // run the specified test, into the tmp dir
  simple_mpi( ref_tree, query_file, qry_info, outdir, options, invocation );

  // read the result
  placement::Sample new_result = placement::JplaceReader().read(
      utils::from_file( outdir + "/epa_result.jplace" ) );
  placement::Sample old_result
      = placement::JplaceReader().read( utils::from_file( previous_jplace ) );

  check_equal( old_result, new_result );
}

TEST( regression, neotrop )
{
  Options options;
  std::string dir = env->data_dir + "neotrop/";

  full_test( dir + "epa_result.jplace",
             dir + "tree.newick",
             dir + "reference.fasta.gz",
             dir + "query_1k.fasta.gz",
             dir + "infofile",
             options );
}