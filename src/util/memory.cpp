#include "util/memory.hpp"

#include "genesis/utils/io/input_stream.hpp"

#include <iomanip>
#include <mutex>

#include "util/maps.hpp"

static size_t lookuptable_footprint( size_t const branches,
                                     size_t const states,
                                     size_t const sites )
{
  size_t size = 0;

  // char_to_posish_ size is static:
  size += 128 * sizeof( size_t );

  // class currently saves mutexes to handle parallel access:
  size += branches * sizeof( std::mutex );

  // the actual lookuptable is, per branch, a matrix of all possible per site likelihoods (pslh) for a given number of states
  // however the number of states includes ambiguity chars:
  auto const effective_states = ( states == 4 ) ? NT_MAP_SIZE : AA_MAP_SIZE;
  auto const pslh_vector      = sites * sizeof( double );

  size += branches * pslh_vector * effective_states;

  return size;
}

static size_t partition_footprint( raxml::Model const& model,
                                   Tree_Numbers const& nums,
                                   int const num_sites,
                                   Options const& options )
{
  size_t size = 0;

  /*
		We make a fake partition, such that we can leave it to pll to figure out cruicial size values like
		states_padded, which we then use in our own calculation.
		Purposefully reusing the function that is used in the Tree constructor here.
	 */

  Tree_Numbers fake_nums;
  fake_nums.tip_nodes   = 1;
  fake_nums.inner_nodes = 1;
  fake_nums.branches    = 1;

  auto partition = make_partition( model,
                                   fake_nums,
                                   num_sites,
                                   options );

  size_t const sites_alloc = partition->asc_additional_sites + partition->sites;

  // eigendecomposition valid
  size += partition->rate_matrices * sizeof( int );

  // number of clv's depends on wether we use tipchars
  size_t num_clvs = nums.inner_nodes * 3;
  if( not( partition->attributes & PLL_ATTRIB_PATTERN_TIP ) ) {
    num_clvs += nums.tip_nodes;
  } else {
    /* reference "create_charmap" in pll.c for this part */
    // account for the charmap
    size += PLL_ASCII_SIZE * sizeof( unsigned char );
    // account fort tipmap
    size += PLL_ASCII_SIZE * sizeof( pll_state_t );

    // ttlookup - upper limit? real calc is pretty difficult, good enough to take the non 4x4 case
    size += 1024 * partition->rate_cats * sizeof( double );

    // account for the tipchars
    size_t const tipchars_buffer = nums.tip_nodes
            * sites_alloc
            * sizeof( unsigned char )
        + nums.tip_nodes * sizeof( unsigned char* ); // account for top level arra
    size += tipchars_buffer;
  }

  size_t const clv_buffer = num_clvs
          * sites_alloc
          * partition->states_padded
          * partition->rate_cats
          * sizeof( double )
      + num_clvs * sizeof( double* ); // account for top level array

  // hypothetical size with log(n) optimization
  auto const log_clv_buffer = log2(num_clvs)
          * sites_alloc
          * partition->states_padded
          * partition->rate_cats
          * sizeof( double )
      + log2(num_clvs) * sizeof( double* ); // account for top level array

  LOG_DBG << "\t\t" << format_byte_num(clv_buffer) << "\tCLV Buffer"
          << " (with log(n) opt: " << format_byte_num(log_clv_buffer) << ")";

  size += clv_buffer;

  // p matrices, they are allocated in a memory saving way (consult pll.c)
  size_t const displacement = ( partition->states_padded - partition->states )
      * ( partition->states_padded )
      * sizeof( double );

  size_t const pmat_buffer = partition->prob_matrices
          * partition->states
          * partition->states_padded
          * partition->rate_cats
          * sizeof( double )
      + displacement
      + partition->prob_matrices * sizeof( double* ); // account for top level array
  size += pmat_buffer;

  // eigenvectors
  size_t const eigen_buffer = partition->rate_matrices
          * partition->states
          * partition->states_padded
          * sizeof( double )
      + partition->rate_matrices * sizeof( double* ); // account for top level array
  size += eigen_buffer;

  // invariant eigenvecs
  size_t const inv_eigen_buffer = partition->rate_matrices
          * partition->states
          * partition->states_padded
          * sizeof( double )
      + partition->rate_matrices * sizeof( double* ); // account for top level array
  size += inv_eigen_buffer;

  // eigenvalues
  size_t const eigenval_buffer = partition->rate_matrices
          * partition->states_padded
          * sizeof( double )
      + partition->rate_matrices * sizeof( double* ); // account for top level array
  size += eigenval_buffer;

  // substitution parameters
  size_t const sub_param_buffer = partition->rate_matrices
          * ( ( partition->states * partition->states - partition->states ) / 2 ) // :thinking:
          * sizeof( double )
      + partition->rate_matrices * sizeof( double* ); // account for top level array
  size += sub_param_buffer;

  // frequencies
  size_t const freq_buffer = partition->rate_matrices
          * partition->states_padded
          * sizeof( double )
      + partition->rate_matrices * sizeof( double* ); // account for top level array
  size += freq_buffer;

  // rates
  size += partition->rate_cats * sizeof( double );
  // rate weights
  size += partition->rate_cats * sizeof( double );
  // proportion of invariant sites
  size += partition->rate_matrices * sizeof( double );
  // site weights
  size += sites_alloc * sizeof( unsigned int );

  // scale buffers
  size_t const scaler_size = ( partition->attributes & PLL_ATTRIB_RATE_SCALERS )
      ? sites_alloc * partition->rate_cats
      : sites_alloc;

  size_t const scaler_buffer = partition->scale_buffers
          * scaler_size
          * sizeof( unsigned int )
      + partition->scale_buffers * sizeof( unsigned int* ); // account for top level array
  size += scaler_buffer;

  pll_partition_destroy( partition );

  return size;
}

size_t estimate_footprint( MSA_Info const& ref_info,
                           MSA_Info const& qry_info,
                           raxml::Model const& model,
                           Options const& options )
{
  if( options.repeats ) {
    LOG_ERR << "Cannot accurately calculate memory footprint when using siterepeats!";
    assert( false );
  }

  size_t size = 0;

  auto const tree_nums = Tree_Numbers( ref_info.sequences() );

  // figure out the true size of the ref alignment
  assert( ref_info.sites() == qry_info.sites() );
  assert( ref_info.gap_mask().size() == qry_info.gap_mask().size() );
  auto const num_sites = options.premasking ? ref_info.nongap_count() : ref_info.sites();

  LOG_DBG << "Memory footprint breakdown:";

  size += partition_footprint( model, tree_nums, num_sites, options );
  LOG_DBG << "\t" << format_byte_num( size ) << "\tPartition Total";

  if( options.prescoring ) {
    auto const lookup_size = lookuptable_footprint( tree_nums.branches, model.num_states(), num_sites );
    LOG_DBG << "\t" << format_byte_num( lookup_size ) << "\tPreplacement Lookup";
    size += lookup_size;
  }

  // account for the overhead induced by the tinytree-encapsulated partitions, which are one per thread


  // size of the fasta input stream buffers
  // currently only one: the query MSA_Stream
  auto const input_stream_buffer_size = genesis::utils::InputStream::BlockLength * 3;
  LOG_DBG << "\t" << format_byte_num( input_stream_buffer_size ) << "\tQuery MSA Inputstream";
  size += input_stream_buffer_size;

  return size;
}

std::string format_byte_num( double size )
{
  constexpr std::array< char const*, 6 > magnitude = { { "", "KiB", "MiB", "GiB", "TiB", "PiB" } };

  size_t lvl = 0;
  while( size > 1024 ) {
    size /= 1024.0;
    lvl++;
  }

  std::ostringstream ret;

  ret << std::setprecision( 1 ) << std::fixed;
  ret << size;
  ret << magnitude[ lvl ];

  return ret.str();
}

std::string format_byte_num( size_t size )
{
  return format_byte_num( static_cast<double>(size) );
}
