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
#include "core/BranchBuffer.hpp"
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
#include "sample/functions.hpp"
#include "seq/MSA.hpp"
#include "tree/Tiny_Tree.hpp"
#include "util/set_manipulators.hpp"
#include "util/Timer.hpp"
#include "util/logging.hpp"
#include "util/stringify.hpp"

#ifdef __MPI
#include "net/epa_mpi_util.hpp"
#endif

using timer = Timer< std::chrono::milliseconds >;

/**
 *  Memsaver mode preplacement.
 *
 * Placement using default pendant lengths, center insertion point
 * (preplacement), done via CLV computation managed by the BranchBuffer.
 */
static void preplace( MSA& msa,
                      BranchBuffer& branchbuf,
                      Sample< Preplacement >& sample,
                      Options const& options,
                      timer* time = nullptr )
{
  LOG_DBG1 << "BranchBuffer based";
#ifdef __OMP
  unsigned int const num_threads
      = options.num_threads ? options.num_threads : omp_get_max_threads();
  omp_set_num_threads( num_threads );
  LOG_DBG1 << "Using threads: " << num_threads;
  LOG_DBG << "Max threads: " << omp_get_max_threads();
#endif

  if( time ) {
    time->start();
  }

  BranchBuffer::container_type branch_chunk;

  while( branchbuf.get_next( branch_chunk ) ) {
    // parallelize over branches: each thread places all queries on its
    // designated branch
#pragma omp parallel for schedule( dynamic )
    for( size_t i = 0; i < branch_chunk.size(); ++i ) {
      auto& branch = branch_chunk[ i ];
      // make shallow copy of tinytree if parallelizing here as well!
      for( size_t seq_id = 0; seq_id < msa.size(); ++seq_id ) {
        sample[ seq_id ][ branch.branch_id() ]
            = branch.preplace( msa[ seq_id ], options );
      }
    }
  }

  if( time ) {
    time->stop();
  }
}

/**
 * Placement using default pendant lengths, center insertion point
 * (preplacement), done via explicit CLV computation (no lookup table).
 */
static void preplace( MSA& msa,
                      Tree& reference_tree,
                      std::vector< pll_unode_t* > const& branches,
                      Sample< Preplacement >& sample,
                      Options const& options,
                      timer* time = nullptr )
{
  if( reference_tree.memsave() ) {
    auto const block_size = options.memory_config.concurrent_branches;
    BranchBuffer branchbuf( &reference_tree, block_size );
    preplace( msa, branchbuf, sample, options, time );
    return;
  }

  LOG_DBG1 << "Full Tree based";

#ifdef __OMP
  unsigned int const num_threads = options.num_threads
      ? options.num_threads
      : omp_get_max_threads();
  omp_set_num_threads( num_threads );
  LOG_DBG1 << "Using threads: " << num_threads;
  LOG_DBG1 << "Max threads: " << omp_get_max_threads();
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

    sample[ seq_id ][ branch_id ] = branch->preplace( msa[ seq_id ], options );

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
static void preplace( MSA& msa,
                      LookupPlacement const& lookup,
                      Sample< Preplacement >& sample,
                      Options const& options,
                      timer* time = nullptr )
{

  LOG_DBG1 << "Lookup based";

#ifdef __OMP
  unsigned int const num_threads
      = options.num_threads ? options.num_threads : omp_get_max_threads();
  omp_set_num_threads( num_threads );
  LOG_DBG1 << "Using threads: " << num_threads;
  LOG_DBG1 << "Max threads: " << omp_get_max_threads();
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

/**
 * BLO placement from a BranchBuffer
 */
static void blo_place( Work const& to_place,
                       MSA& msa,
                       BranchBuffer& branchbuf,
                       Sample< Placement >& sample,
                       Options const& options,
                       size_t const seq_id_offset = 0,
                       timer* time                = nullptr )
{
  LOG_DBG1 << "BranchBuffer based";
#ifdef __OMP
  unsigned int const num_threads
      = options.num_threads ? options.num_threads : omp_get_max_threads();
  omp_set_num_threads( num_threads );
  LOG_DBG1 << "Using threads: " << num_threads;
  LOG_DBG1 << "Max threads: " << omp_get_max_threads();
#else
  unsigned int const num_threads = 1;
#endif

  // work seperately
  if( time ) {
    time->start();
  }

  // split the sample structure such that the parts are thread-local
  std::vector< Sample< Placement > > sample_parts( num_threads );
  // Map from sequence indices to indices in the pquery vector.
  auto seq_lookup_vec
      = std::vector< std::unordered_map< size_t, size_t > >( num_threads );

  // we will want thread-local copies of the tiny-trees that persist for a given
  // thread, if it handles the same branch
  std::vector< std::unique_ptr< Tiny_Tree > > branch_ptrs( num_threads );
  auto prev_branch_id = std::numeric_limits< size_t >::max();

  BranchBuffer::container_type branch_chunk;

  while( branchbuf.get_next( branch_chunk ) ) {
#pragma omp parallel for schedule( dynamic ), firstprivate( prev_branch_id )
    for( size_t i = 0; i < branch_chunk.size(); ++i ) {
      auto& branch         = branch_chunk[ i ];
      auto const branch_id = branch.branch_id();
      // array of all sequence IDs that are to be placed here
      auto const& seqs_of_branch = to_place.at( branch_id );

      // make shallow copy of tinytree if parallelizing here as well!
      for( size_t j = 0; j < seqs_of_branch.size(); ++j ) {
        auto const seq_id = seqs_of_branch[ j ];

#ifdef __OMP
        auto const tid = omp_get_thread_num();
#else
        auto const tid = 0;
#endif
        // get a tiny tree representing the current branch,
        // IF the branch has changed. Overwriting the old variable ensures
        // the now unused previous tiny tree is deallocated
        if( ( branch_id != prev_branch_id ) or not branch_ptrs[ tid ] ) {
          branch_ptrs[ tid ] = std::make_unique< Tiny_Tree >( branch, false );
        }

        auto& local_sample = sample_parts[ tid ];
        auto& seq_lookup   = seq_lookup_vec[ tid ];
        auto const& seq    = msa[ seq_id ];

        if( seq_lookup.count( seq_id ) == 0 ) {
          auto const new_idx
              = local_sample.add_pquery( seq_id_offset + seq_id, seq.header() );
          seq_lookup[ seq_id ] = new_idx;
        }
        assert( seq_lookup.count( seq_id ) > 0 );
        local_sample[ seq_lookup[ seq_id ] ].emplace_back(
            branch_ptrs[ tid ]->blo_place( seq, options ) );
      }
      prev_branch_id = branch_id;
    }
  }

  if( time ) {
    time->stop();
  }
  // merge samples back
  merge( sample, std::move( sample_parts ) );
  collapse( sample );
}

static void blo_place( Work const& to_place,
                       MSA& msa,
                       Tree& reference_tree,
                       std::vector< pll_unode_t* > const& branches,
                       Sample< Placement >& sample,
                       Options const& options,
                       size_t const seq_id_offset = 0,
                       timer* time              = nullptr )
{
  if( reference_tree.memsave() ) {
    auto const block_size = options.memory_config.concurrent_branches;
    BranchBuffer branchbuf( &reference_tree, block_size, to_place );
    blo_place( to_place, msa, branchbuf, sample, options, seq_id_offset, time );
    return;
  }

  LOG_DBG1 << "Full Tree based";

#ifdef __OMP
  unsigned int const num_threads = options.num_threads
      ? options.num_threads
      : omp_get_max_threads();
  omp_set_num_threads( num_threads );
  LOG_DBG1 << "Using threads: " << num_threads;
  LOG_DBG1 << "Max threads: " << omp_get_max_threads();
#else
  unsigned int const num_threads = 1;
#endif

  // split the sample structure such that the parts are thread-local
  std::vector< Sample< Placement > > sample_parts( num_threads );

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
        branch_ptrs[ tid ]->blo_place( seq, options ) );

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

  std::unique_ptr< LookupPlacement > lookup_handler
      = ( options.prescoring and options.memory_config.preplace_lookup_enabled )
      ? std::make_unique< LookupPlacement >( reference_tree, branches, options )
      : nullptr;

  auto reader = make_msa_reader( query_file,
                                 msa_info,
                                 options.premasking,
                                 true );

  size_t num_sequences = 0;
  bool first           = true;
  Work all_work;
  Work blo_work;
  MSA chunk;
  Sample< Preplacement > preplace_result;
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

    if( first or ( num_sequences < options.chunk_size ) ) {

      if( options.prescoring ) {
        preplace_result = Sample< Preplacement >( num_sequences, num_branches );
      } else {
        all_work = Work( std::make_pair( 0, num_branches ),
                         std::make_pair( 0, num_sequences ) );
      }
      first = false;
    }

    if( options.prescoring ) {

      LOG_DBG << "Preplacement" << std::endl;
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

    Sample< Placement > blo_sample;

    LOG_DBG << "BLO Placement" << std::endl;
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
