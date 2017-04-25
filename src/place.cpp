#include "place.hpp"

#include <fstream>
#include <string>
#include <memory>
#include <functional>
#include <limits>

#ifdef __OMP
#include <omp.h>
#endif

#include "file_io.hpp"
#include "jplace_util.hpp"
#include "stringify.hpp"
#include "set_manipulators.hpp"
#include "Log.hpp"
#include "Tiny_Tree.hpp"
#include "mpihead.hpp"
#include "pll_util.hpp"
#include "epa_pll_util.hpp"
#include "Timer.hpp"
#include "Work.hpp"
#include "schedule.hpp"
#include "Lookup_Store.hpp"
#include "Pipeline.hpp"
#include "MSA.hpp"
#include "Work.hpp"
#include "Sample.hpp"

#ifdef __MPI
#include "epa_mpi_util.hpp"
#endif

static std::vector<std::string> split(const std::string & text, const std::string delim)
{
  std::vector<std::string> parts;
  size_t start = 0;
  size_t end = 0;

  do
  {
    end = text.find(delim, start);
    end = (end != std::string::npos) ? end : text.length();
    parts.emplace_back(text.substr(start, end - start));
    start = end + delim.length();
  } while (end != std::string::npos and start <= text.length());

  return parts;
}

static std::string trim(const std::string &s, const char l, const char r)
{
  auto wsfront=std::find_if(s.begin(),s.end(),[l](char c){return c==l;});
  std::advance(wsfront, 1);
  auto wsback=std::find_if(s.rbegin(),s.rend(),[r](char c){return c==r;}).base();
  std::advance(wsback, -1);
  return (wsback<=wsfront ? std::string() : std::string(wsfront,wsback));
}

static void place(const Work& to_place, MSA& msa, Tree& reference_tree,
  const std::vector<pll_utree_t *>& branches, Sample& sample,
  bool do_blo, const Options& options, std::shared_ptr<Lookup_Store>& lookup_store)
{

#ifdef __OMP
  unsigned int num_threads = options.num_threads ? options.num_threads : omp_get_max_threads();
#else
  unsigned int num_threads = 1;
#endif

  // split the sample structure such that the parts are thread-local
  std::vector<Sample> sample_parts(num_threads);
  std::vector<Work> work_parts;
  split(to_place, work_parts, num_threads);

  // work seperately
#ifdef __OMP
  #pragma omp parallel for schedule(dynamic)
#endif
  for (size_t i = 0; i < work_parts.size(); ++i) {
    auto prev_branch_id = std::numeric_limits<size_t>::max();

    std::shared_ptr<Tiny_Tree> branch = nullptr;

    for (auto it : work_parts[i]) {
      auto branch_id = it.branch_id;
      auto seq_id = it.sequence_id;

      if ((branch_id != prev_branch_id) or not branch) {
        branch = std::shared_ptr<Tiny_Tree>(new Tiny_Tree(branches[branch_id], branch_id, reference_tree, do_blo, options, lookup_store));
      }

      sample_parts[i].add_placement(seq_id, branch->place(msa[seq_id]));
      prev_branch_id = branch_id;
    }
  }
  // merge samples back
  merge(sample, sample_parts);
}

static void merge_write_results(const std::string& status_file_name, 
                                const std::string& outdir, 
                                const std::string& newick_string,
                                const std::string& invocation)
{
  std::vector<std::string> part_names;

  // parse back the file with the names, get part names by regex
  std::ifstream status_file(status_file_name);
  std::string line;
  while (std::getline(status_file, line)) {
    merge(part_names, split(trim(line, '[', ']'), ",") );
  }

  for (auto& p : part_names) {
    lgr.dbg() << p << std::endl;
  }

  // create output file
  std::ofstream outfile(outdir + "epa_result.jplace");
  lgr << "\nOutput file: " << outdir + "epa_result.jplace" << std::endl;
  outfile << init_jplace_string(newick_string);
  merge_into(outfile, part_names);
  outfile << finalize_jplace_string(invocation);
  outfile.close();
}

void process( Tree& reference_tree, 
              MSA_Stream& msa_stream, 
              const std::string& outdir,
              const Options& options, 
              const std::string& invocation)
{
  /* ===== COMMON DEFINITIONS ===== */
  int local_rank = 0;

#ifdef __MPI
  int world_size;
  MPI_Comm_rank(MPI_COMM_WORLD, &local_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);

  lgr.dbg() << "World Size: " << world_size << std::endl;

  // shuffle available nodes to the different stages
  schedule_type schedule;
  int local_stage;
  const unsigned int num_stages = options.prescoring ? 4 : 2;

  lgr.dbg() << "Stages: " << num_stages << std::endl;

  unsigned int rebalance = 3;
  unsigned int rebalance_delta = rebalance;
  bool reassign_happened = true;

  std::vector<double> init_diff = options.prescoring ? 
      std::vector<double>{1000.0, 1.0, 1000.0, 1.0} : std::vector<double>{1000.0, 1.0};

  // get initial schedule
  auto init_nps = solve(num_stages, world_size, init_diff);
  assign(local_rank, init_nps, schedule, &local_stage);

  lgr.dbg() << "Schedule: ";
  for (size_t i = 0; i < schedule.size(); ++i) {
    lgr.dbg() << schedule[i].size() << " ";
  }
  lgr.dbg() << std::endl;

  Timer timer;
  Timer dummy;

  const auto EPA_MPI_STAGE_LAST_AGGREGATE
    = options.prescoring ? EPA_MPI_STAGE_2_AGGREGATE : EPA_MPI_STAGE_1_AGGREGATE;
  const auto EPA_MPI_DEDICATED_WRITE_RANK = schedule[EPA_MPI_STAGE_LAST_AGGREGATE][0];

  previous_request_storage_t prev_requests;

#endif // __MPI
  Timer flight_time;
  std::ofstream flight_file(outdir + "stat");

  std::string status_file_name(outdir + "pepa.status");

  lgr << "P-EPA - Massively-Parallel Evolutionary Placement Algorithm" << std::endl;
  lgr << "\nInvocation: \n" << invocation << std::endl;

  const auto chunk_size = options.chunk_size;
  lgr.dbg() << "Chunk size: " << chunk_size << std::endl;

  const auto num_branches = reference_tree.nums().branches;

  // get all edges
  std::vector<pll_utree_t *> branches(num_branches);
  auto num_traversed_branches = utree_query_branches(reference_tree.tree(), &branches[0]);
  if (num_traversed_branches != num_branches) {
    throw std::runtime_error{"Traversing the utree went wrong during pipeline startup!"};
  }

  double lowest = 1e15;
  for (size_t i = 0; i < num_branches; ++i) {
    if(branches[i]->length < lowest)
      lowest = branches[i]->length;
  }

  lgr.dbg() << "smallest BL: " << lowest << std::endl;

  unsigned int chunk_num = 1;
  Sample sample;

  std::shared_ptr<Lookup_Store> previously_calculated_lookups(
    new Lookup_Store(num_branches, reference_tree.partition()->states)
  );

  Work all_work(std::make_pair(0, num_branches), std::make_pair(0, chunk_size));
  Work first_placement_work;
  Work second_placement_work; // dummy structure to be filled during operation
  std::vector<std::string> part_names; // filenames of partial results

#ifdef __MPI
  if (local_rank == 0)
  {
#endif
  std::ofstream trunc_status_file(status_file_name, std::ofstream::trunc);
  flight_time.start();
#ifdef __MPI
  }
#endif
  MSA chunk;
  size_t num_sequences;
  // while data
  // TODO this could be a stream read, such that cat msa.fasta | epa ... would work
  while ((num_sequences = msa_stream.read_next(chunk, chunk_size)) > 0)
  {

    if (num_sequences < chunk_size) {
      all_work = Work(std::make_pair(0, num_branches), std::make_pair(0, num_sequences));
    }

#ifdef __MPI
    // timer.start(); // start timer of any stage
    //==============================================================
    // EPA_MPI_STAGE_1_COMPUTE === BEGIN
    //==============================================================
    if (local_stage == EPA_MPI_STAGE_1_COMPUTE)
    {
    // if previous chunk was a rebalance chunk or this is the first chunk (0 mod anything = 0)
    // ...then we need to correctly assign/reassign the workload of the first compute stage
    if ( reassign_happened )
    {
      lgr.dbg() << "Assigning first stage Work" << std::endl;
      const auto& stage = schedule[EPA_MPI_STAGE_1_COMPUTE];
      std::vector<Work> parts;
      split(all_work, parts, stage.size());
      // find the stage-relative rank
      auto it = std::find(stage.begin(), stage.end(), local_rank);
      size_t stage_rank = std::distance(stage.begin(), it);
      first_placement_work = parts[stage_rank];
      reassign_happened = false;
    }
    timer.start();
#else
    first_placement_work = all_work;
#endif //__MPI

    place(first_placement_work, chunk, reference_tree, branches, sample, !options.prescoring, options, previously_calculated_lookups);

#ifdef __MPI
    timer.stop();
    // MPI: split the result and send the part to correct aggregate node
    lgr.dbg() << "Sending Stage 1 Results..." << std::endl;
    epa_mpi_split_send(sample, schedule[EPA_MPI_STAGE_1_AGGREGATE], MPI_COMM_WORLD, prev_requests, dummy);
    lgr.dbg() << "Stage 1 Send done!" << std::endl;

    } // endif (local_stage == EPA_MPI_STAGE_1_COMPUTE)
    //==============================================================
    // EPA_MPI_STAGE_1_COMPUTE === END
    //==============================================================

    //==============================================================
    // EPA_MPI_STAGE_1_AGGREGATE === BEGIN
    //==============================================================
    if (local_stage == EPA_MPI_STAGE_1_AGGREGATE)
    {
    // (MPI: receive results, merge them)
    lgr.dbg() << "Recieving Stage 1 Results..." << std::endl;
    epa_mpi_receive_merge(sample, schedule[EPA_MPI_STAGE_1_COMPUTE], MPI_COMM_WORLD, dummy);
    lgr.dbg() << "Stage 1 Recieve done!" << std::endl;
    timer.start();
#endif // __MPI

    compute_and_set_lwr(sample);

    // if this was a prescring run, select the candidate edges
    if (options.prescoring) {
      if (options.prescoring_by_percentage) {
        discard_bottom_x_percent(sample, (1.0 - options.prescoring_threshold));
      } else {
        discard_by_accumulated_threshold(sample, options.prescoring_threshold);
      }
    }

    if(options.prescoring)
      second_placement_work = Work(sample);
#ifdef __MPI
    // if prescoring was selected, we need to send the intermediate results off to thorough placement
    if (options.prescoring)
    {
      timer.stop();
      epa_mpi_split_send(second_placement_work, schedule[EPA_MPI_STAGE_2_COMPUTE], MPI_COMM_WORLD, prev_requests, dummy);
    }
    } // endif (local_stage == EPA_MPI_STAGE_1_AGGREGATE)
    //==============================================================
    // EPA_MPI_STAGE_1_AGGREGATE === END
    //==============================================================

    //==============================================================
    // EPA_MPI_STAGE_2_COMPUTE === BEGIN
    //==============================================================
    if (local_stage == EPA_MPI_STAGE_2_COMPUTE and options.prescoring)
    {
    epa_mpi_receive_merge(second_placement_work, schedule[EPA_MPI_STAGE_1_AGGREGATE], MPI_COMM_WORLD, dummy);
    timer.start();
#endif // __MPI
    if (options.prescoring)
    {
      place(second_placement_work, chunk, reference_tree, branches, sample, true, options, previously_calculated_lookups);
    }
#ifdef __MPI
    timer.stop();
    if(options.prescoring) {
      epa_mpi_split_send(sample, schedule[EPA_MPI_STAGE_2_AGGREGATE], MPI_COMM_WORLD, prev_requests, dummy);
    }

    } // endif (local_stage == EPA_MPI_STAGE_2_COMPUTE)
    //==============================================================
    // EPA_MPI_STAGE_2_COMPUTE === END
    //==============================================================

    //==============================================================
    // EPA_MPI_STAGE_2_AGGREGATE === BEGIN
    //==============================================================
    if (local_stage == EPA_MPI_STAGE_LAST_AGGREGATE)
    {
      // only if this is the 4th stage do we need to get from mpi
    if (local_stage == EPA_MPI_STAGE_2_AGGREGATE)
    {
      epa_mpi_receive_merge(sample, schedule[EPA_MPI_STAGE_2_COMPUTE], MPI_COMM_WORLD, dummy);
      timer.start();
    }
#endif // __MPI
    // recompute the lwrs
    if (options.prescoring) {
      compute_and_set_lwr(sample);
    }
    if (options.acc_threshold) {
      lgr.dbg() << "Filtering by accumulated threshold: " << options.support_threshold << std::endl;
      discard_by_accumulated_threshold( sample, 
                                        options.support_threshold,
                                        options.filter_min,
                                        options.filter_max);
    } else {
      lgr.dbg() << "Filtering placements below threshold: " << options.support_threshold << std::endl;
      discard_by_support_threshold( sample,
                                    options.support_threshold,
                                    options.filter_min,
                                    options.filter_max);
    }

    // write results of current last stage aggregator node to a part file
    std::string part_file_name(outdir + "epa." + std::to_string(local_rank)
      + "." + std::to_string(chunk_num) + ".part");
    std::ofstream part_file(part_file_name);
    part_file << sample_to_jplace_string(sample, chunk);
    part_names.push_back(part_file_name);
    part_file.close();

#ifdef __MPI
    // gather file names from all last stage aggregate nodes on foreman, append to status file
    epa_mpi_gather( part_names, 
                    schedule[EPA_MPI_STAGE_LAST_AGGREGATE][0],
                    schedule[EPA_MPI_STAGE_LAST_AGGREGATE],
                    local_rank,
                    dummy
                  );
    if (local_rank == schedule[EPA_MPI_STAGE_LAST_AGGREGATE][0])
    {
#endif
    std::ofstream status_file(status_file_name, std::ofstream::app);
    status_file << chunk_num << ":" << chunk_size << " [";
    for (size_t i = 0; i < part_names.size(); ++i)
    {
      status_file << part_names[i];
      if (i < part_names.size() - 1) status_file << ",";  
    }
    status_file << "]" << std::endl;
    part_names.clear();

#ifdef __MPI
    }
    timer.stop();
    } // endif aggregate cleanup
#endif //__MPI

    lgr.dbg() << "Chunk " << chunk_num << " done!" << std::endl;
    //==============================================================
    // EPA_MPI_STAGE_2_AGGREGATE === END
    //==============================================================
    // timer.stop(); // stop timer of any stage

#ifdef __MPI
    if ( chunk_num == rebalance ) // time to rebalance
    {
      MPI_BARRIER(MPI_COMM_WORLD);
      if (local_rank == 0) 
      {
        flight_time.stop();
        double aft = flight_time.average() / rebalance_delta;
        flight_file << chunk_num << ";" << aft << std::endl;
      }

      lgr.dbg() << "Rebalancing..." << std::endl;
      int foreman = schedule[local_stage][0];
      // Step 0: get per node average
      Timer per_node_avg({timer.average()});
      // Step 1: aggregate the runtime statistics, first at the lowest rank per stage
      lgr.dbg() << "aggregate the runtime statistics..." << std::endl;
      Timer dummy;
      epa_mpi_gather(per_node_avg, foreman, schedule[local_stage], local_rank, dummy);
      lgr.dbg() << "Runtime aggregate done!" << std::endl;

      // Step 2: calculate total time needed per chunk for the stage, reflecting effort spent
      std::vector<double> perstage_total(num_stages);

      int color = (local_rank == foreman) ? 1 : MPI_UNDEFINED;
      MPI_Comm foreman_comm;
      MPI_Comm_split(MPI_COMM_WORLD, color, local_stage, &foreman_comm);

      if (local_rank == foreman)
      {
        auto total_stagetime = per_node_avg.sum();
        // Step 3: make known to all other stage representatives (mpi_allgather)
        lgr.dbg() << "Foremen allgather..." << std::endl;
        MPI_Allgather(&total_stagetime, 1, MPI_DOUBLE, &perstage_total[0], 1, MPI_DOUBLE, foreman_comm);
        lgr.dbg() << "Foremen allgather done!" << std::endl;
        MPI_Comm_free(&foreman_comm);
      }
      // ensure all messages were received and previous requests are cleared
      epa_mpi_waitall(prev_requests);
      
      MPI_BARRIER(MPI_COMM_WORLD);
      // Step 4: stage representatives forward results to all stage members
      // epa_mpi_bcast(perstage_total, foreman, schedule[local_stage], local_rank);
      lgr.dbg() << "Broadcasting..." << std::endl;
      MPI_Comm stage_comm;
      MPI_Comm_split(MPI_COMM_WORLD, local_stage, local_rank, &stage_comm);

      // foreman must always be rank 0 in stage communicator:
      int split_key = local_rank == foreman ? -1 : local_rank;
      MPI_Comm_split(MPI_COMM_WORLD, local_stage, split_key, &stage_comm);
      MPI_Bcast(&perstage_total[0], num_stages, MPI_DOUBLE, 0, stage_comm);
      MPI_Comm_free(&stage_comm);
      lgr.dbg() << "Broadcasting done!" << std::endl;

      lgr.dbg() << "perstage total:";
      for (size_t i = 0; i < perstage_total.size(); ++i)
      {
        lgr.dbg() << " " << perstage_total[i];
      }
      lgr.dbg() << std::endl;

      // Step 5: calculate schedule on every rank, deterministically!
      to_difficulty(perstage_total);

      lgr.dbg() << "perstage difficulty:";
      for (size_t i = 0; i < perstage_total.size(); ++i)
      {
        lgr.dbg() << " " << perstage_total[i];
      }
      lgr.dbg() << std::endl;

      auto nps = solve(num_stages, world_size, perstage_total);
      reassign(local_rank, nps, schedule, &local_stage);
      // Step 6: re-engage pipeline with new assignments
      lgr.dbg() << "New Schedule:";
      for (size_t i = 0; i < nps.size(); ++i)
      {
        lgr.dbg() << " " << nps[i]; 
      }
      lgr.dbg() << std::endl;
      // compute stages should try to keep their edge assignment! affinity!
      lgr.dbg() << "Rebalancing done!" << std::endl;
      // exponential back-off style rebalance:
      rebalance_delta *= 2;
      rebalance += rebalance_delta;
      reassign_happened = true;
      if (local_rank == 0) 
      {
        flight_time.clear();
        flight_time.start();
      }

      prev_requests.clear();
      timer.clear();
    }


#else
    flight_time.stop();
    double aft = flight_time.average();
    flight_file << chunk_num << ";" << aft << std::endl;
    flight_time.clear();
    flight_time.start();
#endif // __MPI
    second_placement_work.clear();
    sample.clear();
    // msa_stream.clear();
    chunk_num++;
  }


  //==============================================================
  // POST COMPUTATION
  //==============================================================
  lgr.dbg() << "Starting Post-Comp" << std::endl;
  // finally, paste all part files together
#ifdef __MPI
  MPI_BARRIER(MPI_COMM_WORLD); // must barrier to avoid paste before finish by rank 0
  if (local_rank == EPA_MPI_DEDICATED_WRITE_RANK)
  {
#endif
    auto newick_string = get_numbered_newick_string(reference_tree.tree());

    merge_write_results(status_file_name, 
                        outdir, 
                        newick_string,
                        invocation);
#ifdef __MPI
  }
#endif
}

void tmp_pipeline_test( Tree& reference_tree, 
                        MSA_Stream& msa_stream, 
                        const std::string& outdir,
                        const Options& options,
                        const std::string& invocation)
{

  auto local_rank = 0;

  lgr << "WARNING! THIS FUNCTION IS EXPERIMENTAL!" << std::endl;

  // Timer flight_time;
  std::ofstream flight_file(outdir + "stat");

  std::string status_file_name(outdir + "pepa.status");
  std::ofstream trunc_status_file(status_file_name, std::ofstream::trunc);

  std::vector<std::string> part_names;

  const auto chunk_size = options.chunk_size;
  lgr.dbg() << "Chunk size: " << chunk_size << std::endl;

  const auto num_branches = reference_tree.nums().branches;

  // get all edges
  std::vector<pll_utree_t *> branches(num_branches);
  auto num_traversed_branches = utree_query_branches(reference_tree.tree(), &branches[0]);
  if (num_traversed_branches != num_branches) {
    throw std::runtime_error{"Traversing the utree went wrong during pipeline startup!"};
  }

  unsigned int chunk_num = 1;
  
  std::shared_ptr<Lookup_Store> lookups(
    new Lookup_Store(num_branches, reference_tree.partition()->states)
  );

  Work all_work(std::make_pair(0, num_branches), std::make_pair(0, chunk_size));
  
  MSA chunk;

  size_t num_sequences = 0;

  auto prehook = [&]() -> void {
    lgr.dbg() << "INGESTING - READING" << std::endl;
    num_sequences = msa_stream.read_next(chunk, chunk_size);
  };

  auto ingestion = [&](VoidToken&) -> Work {
    lgr.dbg() << "INGESTING - CREATING WORK" << std::endl;
    // auto num_sequences = msa_stream.read_next(chunk, chunk_size);
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

  auto preplacement = [&](Work& work) -> Sample {
    lgr.dbg() << "PREPLACING" << std::endl;

    Sample result;

    place(work, 
          chunk, 
          reference_tree, 
          branches, 
          result,
          false, 
          options, 
          lookups
    );
    return result;
  };

  auto candidate_selection = [&](Sample& sample) -> Work {
    lgr.dbg() << "SELECTING CANDIDATES" << std::endl;

    compute_and_set_lwr(sample);

    if (options.prescoring_by_percentage) {
      discard_bottom_x_percent(sample, (1.0 - options.prescoring_threshold));
    } else {
      discard_by_accumulated_threshold(sample, options.prescoring_threshold);
    }
    return Work(sample);
  };


  auto thorough_placement = [&](Work& work) -> Sample {
    lgr.dbg() << "BLO PLACEMENT" << std::endl;

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
    lgr.dbg() << "WRITING" << std::endl;

    compute_and_set_lwr(sample);
    if (options.acc_threshold) {
      lgr.dbg() << "Filtering by accumulated threshold: " << options.support_threshold << std::endl;
      discard_by_accumulated_threshold( sample, 
                                        options.support_threshold,
                                        options.filter_min,
                                        options.filter_max);
    } else {
      lgr.dbg() << "Filtering placements below threshold: " << options.support_threshold << std::endl;
      discard_by_support_threshold( sample,
                                    options.support_threshold,
                                    options.filter_min,
                                    options.filter_max);
    }


    // write results of current last stage aggregator node to a part file
    if (sample.size()) {
      std::string part_file_name(outdir + "epa." + std::to_string(local_rank)
        + "." + std::to_string(chunk_num) + ".part");
      std::ofstream part_file(part_file_name);
      part_file << sample_to_jplace_string(sample, chunk);
      part_names.push_back(part_file_name);
      part_file.close();

      // TODO for MPI, somehow ensure this code is only run on the stage foreman

      std::ofstream status_file(status_file_name, std::ofstream::app);
      status_file << chunk_num << ":" << chunk_size << " [";
      for (size_t i = 0; i < part_names.size(); ++i) {
        status_file << part_names[i];
        if (i < part_names.size() - 1) status_file << ",";  
      }
      status_file << "]" << std::endl;
      part_names.clear();
    }
    
    ++chunk_num;
    return VoidToken();
  };

  auto h_p = make_pipeline(ingestion, prehook)
    .push(preplacement)
    .push(candidate_selection)
    .push(thorough_placement)
    .push(write_result);
  h_p.process();

  auto newick_string = get_numbered_newick_string(reference_tree.tree());

  merge_write_results(status_file_name, 
                      outdir, 
                      newick_string,
                      invocation);

}

