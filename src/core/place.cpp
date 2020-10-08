#include "core/place.hpp"

#include <fstream>
#include <functional>
#include <limits>
#include <memory>
#include <string>

#ifdef __OMP
#include <omp.h>
#endif

#include "core/Lookup_Store.hpp"
#include "core/Work.hpp"
#include "core/heuristics.hpp"
#include "core/pll/epa_pll_util.hpp"
#include "core/pll/pll_util.hpp"
#include "io/Binary_Fasta.hpp"
#include "io/file_io.hpp"
#include "io/jplace_util.hpp"
#include "io/jplace_writer.hpp"
#include "io/msa_reader.hpp"
#include "net/mpihead.hpp"
#include "pipeline/Pipeline.hpp"
#include "pipeline/schedule.hpp"
#include "sample/Sample.hpp"
#include "seq/MSA.hpp"
#include "set_manipulators.hpp"
#include "tree/Tiny_Tree.hpp"
#include "util/Timer.hpp"
#include "util/logging.hpp"
#include "util/stringify.hpp"

#ifdef __MPI
#include "net/epa_mpi_util.hpp"
#endif

using mytimer = Timer< std::chrono::milliseconds >;

/**
 * Placement using default pendant lengths, center insertion point
 * (preplacement), done via explicit CLV computation (no lookup table).
 */
template< class T >
static void preplace( MSA& msa,
                      Tree& reference_tree,
                      std::vector< pll_unode_t* > const& branches,
                      Sample< T >& sample,
                      Options const& options,
                      mytimer* time = nullptr )
{

#ifdef __OMP
  unsigned int const num_threads = options.num_threads
      ? options.num_threads
      : omp_get_max_threads();
  omp_set_num_threads( num_threads );
  LOG_DBG << "Using threads: " << num_threads;
  LOG_DBG << "Max threads: " << omp_get_max_threads();
#else
  unsigned int const num_threads = 1;
#endif

  size_t const num_sequences = msa.size();
  size_t const num_branches  = branches.size();

  std::vector< std::unique_ptr< Tiny_Tree > > branch_ptrs( num_threads );
  auto prev_branch_id = std::numeric_limits< size_t >::max();

  if( time ) {
    time->start();
  }
#ifdef __OMP
#pragma omp parallel for schedule( guided, 10000 ), firstprivate( prev_branch_id )
#endif
  for( size_t i = 0; i < num_sequences * num_branches; ++i ) {

#ifdef __OMP
    auto const tid = omp_get_thread_num();
#else
    auto const tid = 0;
#endif
    // reference to the threadlocal branch
    auto& branch = branch_ptrs[ tid ];

    auto const branch_id = static_cast< size_t >( i ) / num_sequences;
    auto const seq_id    = i % num_sequences;

    // get a tiny tree representing the current branch,
    // IF the branch has changed. Overwriting the old variable ensures
    // the now unused previous tiny tree is deallocated
    if( ( branch_id != prev_branch_id ) or not branch ) {
      // as make_unique produces an rvalue, this is a move assignment and thus legal
      branch = std::make_unique< Tiny_Tree >( branches[ branch_id ],
                                              branch_id,
                                              reference_tree );
    }

    sample[ seq_id ][ branch_id ] = branch->place(
        msa[ seq_id ],
        false,
        options );

    prev_branch_id = branch_id;
  }
  if( time ) {
    time->stop();
  }
}

/**
 * Lookup-based placement using default pendant lengths, center insertion point
 * (preplacement).
 */
template< class T >
static void preplace( MSA& msa,
                      LookupPlacement const& lookup,
                      Sample< T >& sample,
                      Options const& options,
                      mytimer* time = nullptr )
{

#ifdef __OMP
  unsigned int const num_threads
      = options.num_threads ? options.num_threads : omp_get_max_threads();
  omp_set_num_threads( num_threads );
  LOG_DBG << "Using threads: " << num_threads;
  LOG_DBG << "Max threads: " << omp_get_max_threads();
#endif

  auto const num_sequences = msa.size();
  auto const num_branches  = lookup.num_branches();

  if( time ) {
    time->start();
  }

// TODO: currently uses the old way of accessing "Work" packages, double for
// loop would be easier to undersand, possibly more efficient by OMP
#ifdef __OMP
#pragma omp parallel for schedule( guided, 10000 )
#endif
  for( size_t i = 0; i < num_sequences * num_branches; ++i ) {

    auto const branch_id = static_cast< size_t >( i ) / num_sequences;
    auto const seq_id    = i % num_sequences;

    sample[ seq_id ][ branch_id ]
        = lookup.place( branch_id, msa[ seq_id ], options.premasking );
  }
  if( time ) {
    time->stop();
  }
}

template< class T >
static void blo_place( Work const& to_place,
                       MSA& msa,
                       Tree& reference_tree,
                       std::vector< pll_unode_t* > const& branches,
                       Sample< T >& sample,
                       Options const& options,
                       size_t const seq_id_offset = 0,
                       mytimer* time              = nullptr )
{

#ifdef __OMP
  unsigned int const num_threads = options.num_threads
      ? options.num_threads
      : omp_get_max_threads();
  omp_set_num_threads( num_threads );
  LOG_DBG << "Using threads: " << num_threads;
  LOG_DBG << "Max threads: " << omp_get_max_threads();
#else
  unsigned int const num_threads = 1;
#endif

  // split the sample structure such that the parts are thread-local
  std::vector< Sample< T > > sample_parts( num_threads );

  // build vector of elements
  std::vector< Work::Work_Pair > id;
  for( auto it = to_place.begin(); it != to_place.end(); ++it ) {
    id.push_back( *it );
  }

  // Map from sequence indices to indices in the pquery vector.
  auto seq_lookup_vec = std::vector< std::unordered_map< size_t, size_t > >( num_threads );

  std::vector< std::unique_ptr< Tiny_Tree > > branch_ptrs( num_threads );
  auto prev_branch_id = std::numeric_limits< size_t >::max();

  // work seperately
  if( time ) {
    time->start();
  }
#ifdef __OMP
#pragma omp parallel for schedule( dynamic ), firstprivate( prev_branch_id )
#endif
  for( size_t i = 0; i < id.size(); ++i ) {

#ifdef __OMP
    auto const tid = omp_get_thread_num();
#else
    auto const tid               = 0;
#endif
    auto& local_sample = sample_parts[ tid ];
    auto& seq_lookup   = seq_lookup_vec[ tid ];

    auto const branch_id = id[ i ].branch_id;
    auto const seq_id    = id[ i ].sequence_id;
    auto const& seq      = msa[ seq_id ];

    // get a tiny tree representing the current branch,
    // IF the branch has changed. Overwriting the old variable ensures
    // the now unused previous tiny tree is deallocated
    if( ( branch_id != prev_branch_id ) or not branch_ptrs[ tid ] ) {
      // as make_unique produces an rvalue, this is a move assignment and thus legal
      branch_ptrs[ tid ] = std::make_unique< Tiny_Tree >( branches[ branch_id ],
                                                          branch_id,
                                                          reference_tree );
    }

    if( seq_lookup.count( seq_id ) == 0 ) {
      auto const new_idx   = local_sample.add_pquery( seq_id_offset + seq_id, seq.header() );
      seq_lookup[ seq_id ] = new_idx;
    }
    assert( seq_lookup.count( seq_id ) > 0 );
    local_sample[ seq_lookup[ seq_id ] ].emplace_back(
        branch_ptrs[ tid ]->place( seq,
                                   true,
                                   options ) );

    prev_branch_id = branch_id;
  }
  if( time ) {
    time->stop();
  }
  // merge samples back
  merge( sample, std::move( sample_parts ) );
  collapse( sample );
}

void simple_mpi( Tree& reference_tree,
                 std::string const& query_file,
                 MSA_Info const& msa_info,
                 std::string const& outdir,
                 Options const& options,
                 std::string const& invocation )
{
  auto const num_branches = reference_tree.nums().branches;

  // get all edges
  std::vector< pll_unode_t* > branches( num_branches );
  auto num_traversed_branches = utree_query_branches( reference_tree.tree(), &branches[ 0 ] );
  if( num_traversed_branches != num_branches ) {
    throw std::runtime_error{ "Traversing the utree went wrong during pipeline startup!" };
  }

  std::unique_ptr< LookupPlacement > lookup_handler = ( options.preplacement_lookup )
      ? std::make_unique< LookupPlacement >( reference_tree, branches )
      : nullptr;

  auto reader = make_msa_reader( query_file,
                                 msa_info,
                                 options.premasking,
                                 true );

  using Sample         = Sample< Placement >;
  size_t num_sequences = 0;
  bool first           = true;
  Work all_work;
  Work blo_work;
  MSA chunk;
  Sample preplace_result;
  size_t sequences_done = 0; // not just for info output!

  // prepare output file
  LOG_INFO << "Output file: " << outdir + "epa_result.jplace";
  jplace_writer jplace( outdir, "epa_result.jplace",
                        get_numbered_newick_string( reference_tree.tree(),
                                                    reference_tree.mapper(),
                                                    options.precision ),
                        invocation,
                        reference_tree.mapper() );
  jplace.set_precision( options.precision );

  while( ( num_sequences = reader->read_next( chunk, options.chunk_size ) ) ) {

    assert( chunk.size() == num_sequences );

    LOG_DBG << "num_sequences: " << num_sequences << std::endl;

    size_t const seq_id_offset = sequences_done + reader->local_seq_offset();

    if( first or num_sequences < options.chunk_size ) {
      all_work = Work( std::make_pair( 0, num_branches ), std::make_pair( 0, num_sequences ) );

      if( options.prescoring ) {
        preplace_result = Sample( num_sequences, num_branches );
      }
      first = false;
    }

    if( options.prescoring ) {

      LOG_DBG << "Preplacement." << std::endl;
      if( lookup_handler ) {
        // if we have a preplacement lookup, use that
        preplace( chunk,
                  *lookup_handler,
                  preplace_result,
                  options );
      } else {
        // otherwise do placement the hard way: through CLV calculation
        preplace( chunk,
                  reference_tree,
                  branches,
                  preplace_result,
                  options );
      }

      LOG_DBG << "Selecting candidates." << std::endl;

      blo_work = apply_heuristic( preplace_result, options );

    } else {
      blo_work = all_work;
    }

    Sample blo_sample;

    LOG_DBG << "BLO Placement." << std::endl;
    blo_place( blo_work,
               chunk,
               reference_tree,
               branches,
               blo_sample,
               options,
               seq_id_offset );

    // Output
    compute_and_set_lwr( blo_sample );
    filter( blo_sample, options );

    // pass the result chunk to the writer
    jplace.write( blo_sample );

    sequences_done += num_sequences;
    LOG_INFO << sequences_done << " Sequences done!";
  }

  jplace.wait();

  MPI_BARRIER( MPI_COMM_WORLD );
}
