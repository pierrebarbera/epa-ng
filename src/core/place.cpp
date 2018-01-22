#include "core/place.hpp"

#include <fstream>
#include <string>
#include <memory>
#include <functional>
#include <limits>

#ifdef __OMP
#include <omp.h>
#endif

#include "io/file_io.hpp"
#include "io/jplace_util.hpp"
#include "io/msa_reader.hpp"
#include "io/Binary_Fasta.hpp"
#include "io/jplace_writer.hpp"
#include "util/stringify.hpp"
#include "util/logging.hpp"
#include "util/Timer.hpp"
#include "tree/Tiny_Tree.hpp"
#include "net/mpihead.hpp"
#include "pipeline/schedule.hpp"
#include "pipeline/Pipeline.hpp"
#include "seq/MSA.hpp"
#include "core/pll/pll_util.hpp"
#include "core/pll/epa_pll_util.hpp"
#include "core/Work.hpp"
#include "core/Lookup_Store.hpp"
#include "core/Work.hpp"
#include "sample/Sample.hpp"
#include "set_manipulators.hpp"

#ifdef __MPI
#include "net/epa_mpi_util.hpp"
#endif

template <class T>
static void place(const Work& to_place,
                  MSA& msa,
                  Tree& reference_tree,
                  const std::vector<pll_unode_t *>& branches,
                  Sample<T>& sample,
                  bool do_blo,
                  const Options& options,
                  std::shared_ptr<Lookup_Store>& lookup_store,
                  const size_t seq_id_offset=0,
                  Timer<>* time=nullptr)
{

#ifdef __OMP
  const unsigned int num_threads  = options.num_threads
                                  ? options.num_threads
                                  : omp_get_max_threads();
  omp_set_num_threads(num_threads);
  LOG_DBG << "Using threads: " << num_threads;
  LOG_DBG << "Max threads: " << omp_get_max_threads();
#else
  const unsigned int num_threads = 1;
#endif

  // split the sample structure such that the parts are thread-local
  std::vector<Sample<T>> sample_parts(num_threads);

  // build vector of elements
  std::vector<Work::Work_Pair> id;
  for(auto it = to_place.begin(); it != to_place.end(); ++it) {
    id.push_back(*it);
  }

  // Map from sequence indices to indices in the pquery vector.
  auto seq_lookup_vec = std::vector<std::unordered_map<size_t, size_t>>(num_threads);

  std::vector<std::unique_ptr<Tiny_Tree>> branch_ptrs(num_threads);
  auto prev_branch_id = std::numeric_limits<size_t>::max();
  
  // work seperately
  if (time){
    time->start();
  }
#ifdef __OMP
  #pragma omp parallel for schedule(guided, 5000), firstprivate(prev_branch_id)
#endif
  for (size_t i = 0; i < id.size(); ++i) {

#ifdef __OMP
    const auto tid = omp_get_thread_num();
#else
    const auto tid = 0;
#endif
    auto& local_sample = sample_parts[tid];
    auto& seq_lookup = seq_lookup_vec[tid];

    const auto branch_id = id[i].branch_id;
    const auto seq_id = id[i].sequence_id;
    const auto& seq = msa[seq_id];

    // get a tiny tree representing the current branch,
    // IF the branch has changed. Overwriting the old variable ensures
    // the now unused previous tiny tree is deallocated
    if ((branch_id != prev_branch_id) or not branch_ptrs[tid]) {
      // as make_unique produces an rvalue, this is a move assignment and thus legal
      branch_ptrs[tid] = std::make_unique<Tiny_Tree>(branches[branch_id],
                                           branch_id,
                                           reference_tree,
                                           do_blo,
                                           options,
                                           lookup_store);
    }

    if (seq_lookup.count( seq_id ) == 0) {
      auto const new_idx = local_sample.add_pquery( seq_id_offset + seq_id, seq.header() );
      seq_lookup[ seq_id ] = new_idx;
    }
    assert( seq_lookup.count( seq_id ) > 0 );
    local_sample[ seq_lookup[ seq_id ] ].emplace_back( branch_ptrs[tid]->place(seq) );

    prev_branch_id = branch_id;
  }
  if (time){
    time->stop();
  }
  // merge samples back
  merge(sample, std::move(sample_parts));
  collapse(sample);
}

template <class T>
static void place_thorough(const Work& to_place,
                  MSA& msa,
                  Tree& reference_tree,
                  const std::vector<pll_unode_t *>& branches,
                  Sample<T>& sample,
                  bool,
                  const Options& options,
                  std::shared_ptr<Lookup_Store>& lookup_store,
                  const size_t seq_id_offset=0,
                  Timer<>* time=nullptr)
{
  
#ifdef __OMP
  const unsigned int num_threads  = options.num_threads
                                  ? options.num_threads
                                  : omp_get_max_threads();
  omp_set_num_threads(num_threads);
  LOG_DBG << "Using threads: " << num_threads;
  LOG_DBG << "Max threads: " << omp_get_max_threads();
#else
  const unsigned int num_threads = 1;
#endif

  // split the sample structure such that the parts are thread-local
  std::vector<Sample<T>> sample_parts(num_threads);

  // build vector of elements
  std::vector<Work::Work_Pair> id;
  for(auto it = to_place.begin(); it != to_place.end(); ++it) {
    id.push_back(*it);
  }

  // Map from sequence indices to indices in the pquery vector.
  auto seq_lookup_vec = std::vector<std::unordered_map<size_t, size_t>>(num_threads);

  std::vector<std::unique_ptr<Tiny_Tree>> branch_ptrs(num_threads);
  auto prev_branch_id = std::numeric_limits<size_t>::max();
  
  // work seperately
  if (time){
    time->start();
  }
#ifdef __OMP
  #pragma omp parallel for schedule(dynamic), firstprivate(prev_branch_id)
#endif
  for (size_t i = 0; i < id.size(); ++i) {

#ifdef __OMP
    const auto tid = omp_get_thread_num();
#else
    const auto tid = 0;
#endif
    auto& local_sample = sample_parts[tid];
    auto& seq_lookup = seq_lookup_vec[tid];

    const auto branch_id = id[i].branch_id;
    const auto seq_id = id[i].sequence_id;
    const auto& seq = msa[seq_id];

    // get a tiny tree representing the current branch,
    // IF the branch has changed. Overwriting the old variable ensures
    // the now unused previous tiny tree is deallocated
    if ((branch_id != prev_branch_id) or not branch_ptrs[tid]) {
      // as make_unique produces an rvalue, this is a move assignment and thus legal
      branch_ptrs[tid] = std::make_unique<Tiny_Tree>(branches[branch_id],
                                           branch_id,
                                           reference_tree,
                                           true,
                                           options,
                                           lookup_store);
    }

    if (seq_lookup.count( seq_id ) == 0) {
      auto const new_idx = local_sample.add_pquery( seq_id_offset + seq_id, seq.header() );
      seq_lookup[ seq_id ] = new_idx;
    }
    assert( seq_lookup.count( seq_id ) > 0 );
    local_sample[ seq_lookup[ seq_id ] ].emplace_back( branch_ptrs[tid]->place(seq) );

    prev_branch_id = branch_id;
  }
  if (time){
    time->stop();
  }
  // merge samples back
  merge(sample, std::move(sample_parts));
  collapse(sample);
}

void pipeline_place(Tree& reference_tree,
                    const std::string& query_file,
                    const std::string& outdir,
                    const Options& options,
                    const std::string& invocation)
{
  LOG_INFO << "WARNING! THIS FUNCTION IS EXPERIMENTAL!" << std::endl;

  // Timer<> flight_time;
  std::ofstream flight_file(outdir + "stat");

  std::string status_file_name(outdir + "pepa.status");
  std::ofstream trunc_status_file(status_file_name, std::ofstream::trunc);

  std::vector<std::string> part_names;

  const auto chunk_size = options.chunk_size;
  LOG_DBG << "Chunk size: " << chunk_size << std::endl;

  const auto num_branches = reference_tree.nums().branches;

  // get all edges
  std::vector<pll_unode_t *> branches(num_branches);
  auto num_traversed_branches = utree_query_branches(reference_tree.tree(), &branches[0]);
  if (num_traversed_branches != num_branches) {
    throw std::runtime_error{"Traversing the utree went wrong during pipeline startup!"};
  }

  unsigned int chunk_num = 0;
  
  auto lookups = 
    std::make_shared<Lookup_Store>(num_branches, reference_tree.partition()->states);

  Work all_work(std::make_pair(0, num_branches), std::make_pair(0, chunk_size));
  
  MSA chunk;
  Binary_Fasta_Reader reader(query_file);

  size_t num_sequences = 0;

  // create output file
  std::ofstream outfile;

  using Slim_Sample = Sample<Slim_Placement>;
  using Sample      = Sample<Placement>;
  
  // ============ LAMBDAS ============================
  
  // only on one rank, only once at the beginning of pipeline
  auto init_pipe_func = [&]() -> void {
    outfile.open(outdir + "epa_result.jplace");
    auto newick_string = get_numbered_newick_string(reference_tree.tree());
    outfile << init_jplace_string(newick_string);
  };

  auto perloop_prehook = [&]() -> void {
    LOG_DBG << "INGESTING - READING" << std::endl;
    num_sequences = reader.read_next(chunk, chunk_size);
    ++chunk_num;
  };

  auto ingestion = [&](VoidToken&) -> Work {
    LOG_DBG << "INGESTING - CREATING WORK" << std::endl;
    if (num_sequences <= 0) {
      Work work;
      work.is_last(true);
      return work;
    } else if (num_sequences < chunk_size) {
      return Work(std::make_pair(0, num_branches), std::make_pair(0, num_sequences));
    } else {
      return all_work;
    }
  };

  auto preplacement = [&](Work& work) -> Slim_Sample {
    LOG_DBG << "PREPLACING" << std::endl;

    Slim_Sample result;

    place(work,
          chunk,
          reference_tree,
          branches,
          result,
          false,
          options,
          lookups);

    return result;
  };

  // auto ingest_preplace = [&](VoidToken&) -> Slim_Sample {
  //   LOG_DBG << "INGEST & PREPLACE" << std::endl;

  //   Slim_Sample result;

  //   if (num_sequences <= 0) {
  //     result.is_last(true);
  //     return result;
  //   } else if (num_sequences < chunk_size) {
  //     Work work(std::make_pair(0, num_branches), std::make_pair(0, num_sequences));
  //     std::vector<Work> parts
  //   } else {
  //     return all_work;
  //   }


  //   place(work,
  //         chunk,
  //         reference_tree,
  //         branches,
  //         result,
  //         false,
  //         options,
  //         lookups);

  //   return result;
  // };

  auto candidate_selection = [&](Slim_Sample& slim) -> Work {
    LOG_DBG << "SELECTING CANDIDATES" << std::endl;

    Sample sample(slim);

    compute_and_set_lwr(sample);

    if (options.prescoring_by_percentage) {
      discard_bottom_x_percent(sample, 
                              (1.0 - options.prescoring_threshold));
    } else {
      discard_by_accumulated_threshold(sample, 
                                      options.prescoring_threshold,
                                      options.filter_min,
                                      options.filter_max);
    }

    return Work(sample);
  };


  auto thorough_placement = [&](Work& work) -> Sample {
    LOG_DBG << "BLO PLACEMENT" << std::endl;

    Sample result;
    place(work,
          chunk,
          reference_tree,
          branches,
          result,
          true,
          options,
          lookups
    );
    return result;
  };

  auto write_result = [&](Sample& sample) -> VoidToken {
    LOG_DBG << "WRITING" << std::endl;

    compute_and_set_lwr(sample);
    filter(sample, options);

    // write results of current last stage aggregator node to a part file
    if (sample.size()) {
      if (chunk_num > 1) {
        outfile << ",";
      }
      outfile << sample_to_jplace_string(sample);
      // std::string part_file_name(outdir + "epa." + std::to_string(local_rank)
      //   + "." + std::to_string(chunk_num) + ".part");
      // std::ofstream part_file(part_file_name);
      // part_file << sample_to_jplace_string(sample, chunk);
      // part_names.push_back(part_file_name);
      // part_file.close();

      // // TODO for MPI, somehow ensure this code is only run on the stage foreman

      // std::ofstream status_file(status_file_name, std::ofstream::app);
      // status_file << chunk_num << ":" << chunk_size << " [";
      // for (size_t i = 0; i < part_names.size(); ++i) {
      //   status_file << part_names[i];
      //   if (i < part_names.size() - 1) status_file << ",";  
      // }
      // status_file << "]" << std::endl;
      // part_names.clear();
    }

    LOG_INFO << chunk_num * chunk_size  << " Sequences done!"; 

    return VoidToken();
  };

  // only on one rank, only once at the end of the pipeline
  auto finalize_pipe_func = [&]() -> void {
    LOG_INFO << "Output file: " << outdir + "epa_result.jplace";
    outfile << finalize_jplace_string(invocation);
    outfile.close();
  };



  if (options.prescoring) {
    auto pipe = make_pipeline(ingestion, perloop_prehook, init_pipe_func, finalize_pipe_func)
      .push(preplacement)
      .push(candidate_selection)
      .push(thorough_placement)
      .push(write_result);

    pipe.process();
  } else {
    auto pipe = make_pipeline(ingestion, perloop_prehook, init_pipe_func, finalize_pipe_func)
      .push(thorough_placement)
      .push(write_result);

    pipe.process();
  }
}

void simple_mpi(Tree& reference_tree,
                const std::string& query_file,
                const MSA_Info& msa_info,
                const std::string& outdir,
                const Options& options,
                const std::string& invocation)
{
  // Timer<> preplacement_timer;
  // Timer<> preplacement_core_timer;
  // Timer<> thorough_timer;
  // Timer<> thorough_core_timer;

  const auto num_branches = reference_tree.nums().branches;

  // get all edges
  std::vector<pll_unode_t *> branches(num_branches);
  auto num_traversed_branches = utree_query_branches(reference_tree.tree(), &branches[0]);
  if (num_traversed_branches != num_branches) {
    throw std::runtime_error{"Traversing the utree went wrong during pipeline startup!"};
  }

  auto lookups = 
    std::make_shared<Lookup_Store>(num_branches, reference_tree.partition()->states);

  // some MPI prep
  int local_rank = 0;
  int num_ranks = 1;

  MPI_COMM_RANK(MPI_COMM_WORLD, &local_rank);
  MPI_COMM_SIZE(MPI_COMM_WORLD, &num_ranks);

  LOG_INFO << "Number of ranks: " << num_ranks;

  std::vector<int> all_ranks(num_ranks);
  for (int i = 0; i < num_ranks; ++i) {
    all_ranks[i] = i;
  }

  auto reader = make_msa_reader(query_file, msa_info, options.premasking);

  size_t local_rank_seq_offset = 0;

  if (num_ranks > 1) {
    // how many should each rank read?
    const size_t part_size = ceil(reader->num_sequences() / static_cast<double>(num_ranks));
    LOG_INFO << "Number of sequences per rank: " << part_size;

    // read only the locally relevant part of the queries
    // ... by skipping the appropriate amount
    local_rank_seq_offset = part_size * local_rank;
    reader->skip_to_sequence( local_rank_seq_offset );
    // and limiting the reading to the given window
    reader->constrain(part_size);
  }

  size_t num_sequences = 0;
  Work all_work(std::make_pair(0, num_branches), std::make_pair(0, options.chunk_size));

  Work blo_work;

  size_t chunk_num = 1;

  using Sample = Sample<Placement>;
  std::future<void> prev_gather;
  MSA chunk;
  size_t sequences_done = 0; // not just for info output!

  // prepare output file
  Jplace_writer jplace;

  if (local_rank == 0) {
    LOG_INFO << "Output file: " << outdir + "epa_result.jplace";
    jplace = Jplace_writer( outdir + "epa_result.jplace",
                            get_numbered_newick_string(reference_tree.tree()),
                            invocation);
  }

  while ( (num_sequences = reader->read_next(chunk, options.chunk_size)) ) {

    assert(chunk.size() == num_sequences);

    LOG_DBG << "num_sequences: " << num_sequences << std::endl;

    const size_t seq_id_offset = sequences_done + local_rank_seq_offset;

    if (num_sequences < options.chunk_size) {
      all_work = Work(std::make_pair(0, num_branches), std::make_pair(0, num_sequences));
    }

    if (options.prescoring) {

      Sample preplace;

      LOG_DBG << "Preplacement." << std::endl;
      place(all_work,
            chunk,
            reference_tree,
            branches,
            preplace,
            false,
            options,
            lookups,0);

      // Candidate Selection
      LOG_DBG << "Selecting candidates." << std::endl;
      compute_and_set_lwr(preplace);

      if (options.prescoring_by_percentage) {
        discard_bottom_x_percent(preplace, 
                                (1.0 - options.prescoring_threshold));
      } else {
        discard_by_accumulated_threshold( preplace, 
                                          options.prescoring_threshold,
                                          options.filter_min,
                                          options.filter_max);
      }

      blo_work = Work(preplace);

    } else {
      blo_work = all_work;
    }

    Sample blo_sample;

    // BLO placement
    LOG_DBG << "BLO Placement." << std::endl;
    place_thorough( blo_work,
                    chunk,
                    reference_tree,
                    branches,
                    blo_sample,
                    true,
                    options,
                    lookups,
                    seq_id_offset);

    // Output
    compute_and_set_lwr(blo_sample);
    filter(blo_sample, options);

    // pass the result chunk to the writer
    jplace.gather_write(blo_sample, all_ranks, local_rank);

    sequences_done += num_sequences;
    LOG_INFO << sequences_done  << " Sequences done!";
    ++chunk_num;
  }

  jplace.wait();

  MPI_BARRIER(MPI_COMM_WORLD);
}

