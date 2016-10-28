#include "place.hpp"

#include <fstream>
#include <chrono>
#include <omp.h>

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

#ifdef __MPI
#include "epa_mpi_util.hpp"
#endif

static void place(const Work& to_place, MSA_Stream& msa, Tree& reference_tree, 
  const std::vector<pll_utree_t *>& branches, Sample& sample, bool thorough)
{

#ifdef __OMP
  unsigned int num_threads = omp_get_max_threads();
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
  for (size_t i = 0; i < work_parts.size(); ++i)
  {
    for (const auto& pair : work_parts[i])
    {
      auto branch_id = pair.first;
      auto branch = Tiny_Tree(branches[branch_id], branch_id, reference_tree, thorough);
      for (const auto& seq_id : pair.second) 
        sample_parts[i].add_placement(seq_id, branch.place(msa[seq_id]));
    }
  }
  // merge samples back
  merge(sample, sample_parts);
}

void process(Tree& reference_tree, MSA_Stream& msa_stream, const std::string& outdir,
              const Options& options, const std::string& invocation)
{
  /* ===== COMMON DEFINITIONS ===== */
  int local_rank = 0;

#ifdef __MPI
  int world_size;
  MPI_Comm_rank(MPI_COMM_WORLD, &local_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);

  lgr.dbg() << "World Size: " << world_size << std::endl;

  // shuffle available nodes to the different stages
  std::vector<std::vector<int>> schedule;
  int local_stage;
  const unsigned int num_stages = options.prescoring ? 4 : 2;

  lgr.dbg() << "Stages: " << num_stages << std::endl;

  unsigned int rebalance = 10;

  std::vector<double> init_diff = options.prescoring 
    ? std::vector<double>{1000.0, 1.0, 1000.0, 1.0} : std::vector<double>{1000.0, 1.0};

  // get initial schedule
  auto init_nps = solve(num_stages, world_size, init_diff);
  assign(local_rank, init_nps, schedule, &local_stage);

  lgr.dbg() << "Schedule: ";
  for (size_t i = 0; i < schedule.size(); ++i)
  {
    lgr.dbg() << schedule[i].size() << " ";
  }
  lgr.dbg() << std::endl;

  Timer timer;

  const auto EPA_MPI_STAGE_LAST_AGGREGATE 
    = options.prescoring ? EPA_MPI_STAGE_2_AGGREGATE : EPA_MPI_STAGE_1_AGGREGATE;
  const auto EPA_MPI_DEDICATED_WRITE_RANK = schedule[EPA_MPI_STAGE_LAST_AGGREGATE][0];

#endif // __MPI
  lgr << "EPA - Evolutionary Placement Algorithm\n";
  lgr << "\nInvocation: \n" << invocation << "\n";

  const unsigned int chunk_size = 100;
  unsigned int num_sequences;

  const auto num_branches = reference_tree.nums().branches;

  // get all edges
  std::vector<pll_utree_t *> branches(num_branches);
  auto num_traversed_branches = utree_query_branches(reference_tree.tree(), &branches[0]);
  assert(num_traversed_branches == num_branches);

  unsigned int chunk_num = 1;
  Sample sample;

  Work all_work(std::make_pair(0, num_branches), std::make_pair(0, chunk_size));
  Work first_placement_work;
  Work second_placement_work; // dummy structure to be filled during operation
  std::vector<std::string> part_names; // filenames of partial results

  // while data
  // TODO this could be a stream read, such that cat msa.fasta | epa ... would work
  while ((num_sequences = msa_stream.read_next(chunk_size)) > 0)
  {
    if (num_sequences < chunk_size)
      all_work = Work(std::make_pair(0, num_branches), std::make_pair(0, num_sequences));

#ifdef __MPI
    timer.start(); // start timer of any stage
    //==============================================================
    // EPA_MPI_STAGE_1_COMPUTE === BEGIN
    //==============================================================
    if (local_stage == EPA_MPI_STAGE_1_COMPUTE)
    {

    // if previous chunk was a rebalance chunk or this is the first chunk (0 mod anything = 0)
    // ...then we need to correctly assign/reassign the workload of the first compute stage
    if ( !((chunk_num - 1) % rebalance) )
    { 
      lgr.dbg() << "Assigning first stage Work" << std::endl;
      const auto& stage = schedule[EPA_MPI_STAGE_1_COMPUTE];
      std::vector<Work> parts;
      split(all_work, parts, stage.size());
      // find the stage-relative rank
      auto it = std::find(stage.begin(), stage.end(), local_rank);
      size_t stage_rank = std::distance(stage.begin(), it);
      first_placement_work = parts[stage_rank];
    }
#else
    first_placement_work = all_work;
#endif //__MPI    
    
    place(first_placement_work, msa_stream, reference_tree, branches, sample, !options.prescoring);

#ifdef __MPI
    // MPI: split the result and send the part to correct aggregate node
    lgr.dbg() << "Sending Stage 1 Results..." << std::endl;
    epa_mpi_split_send(sample, schedule[EPA_MPI_STAGE_1_AGGREGATE], MPI_COMM_WORLD);
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
    // (MPI: recieve results, merge them)
    lgr.dbg() << "Recieving Stage 1 Results..." << std::endl;
    epa_mpi_recieve_merge(sample, schedule[EPA_MPI_STAGE_1_COMPUTE], MPI_COMM_WORLD);
    lgr.dbg() << "Stage 1 Recieve done!" << std::endl;
#endif // __MPI

    compute_and_set_lwr(sample);

    // if this was a prescring run, select the candidate edges
    if (options.prescoring)
    {
      if (options.prescoring_by_percentage)
        discard_bottom_x_percent(sample, (1.0 - options.prescoring_threshold));
      else
        discard_by_accumulated_threshold(sample, options.prescoring_threshold);
    }
    
    if(options.prescoring)
      second_placement_work = Work(sample);
#ifdef __MPI
    // if prescoring was selected, we need to send the intermediate results off to thorough placement
    if (options.prescoring)
      epa_mpi_split_send(second_placement_work, schedule[EPA_MPI_STAGE_2_COMPUTE], MPI_COMM_WORLD);
    } // endif (local_stage == EPA_MPI_STAGE_1_AGGREGATE)
    //==============================================================
    // EPA_MPI_STAGE_1_AGGREGATE === END
    //==============================================================

    //==============================================================
    // EPA_MPI_STAGE_2_COMPUTE === BEGIN
    //==============================================================
    if (local_stage == EPA_MPI_STAGE_2_COMPUTE and options.prescoring)
    {
    epa_mpi_recieve_merge(second_placement_work, schedule[EPA_MPI_STAGE_1_AGGREGATE], MPI_COMM_WORLD);
#endif // __MPI
    if (options.prescoring)
    {
      place(second_placement_work, msa_stream, reference_tree, branches, sample, true);
    }
#ifdef __MPI
    if(options.prescoring)
      epa_mpi_split_send(sample, schedule[EPA_MPI_STAGE_2_AGGREGATE], MPI_COMM_WORLD);

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
        epa_mpi_recieve_merge(sample, schedule[EPA_MPI_STAGE_2_COMPUTE], MPI_COMM_WORLD);

#endif // __MPI
    // recompute the lwrs
    if (options.prescoring)
      compute_and_set_lwr(sample);
    // discard uninteresting placements
    if (options.acc_threshold)
      discard_by_accumulated_threshold(sample, options.support_threshold);
    else
      discard_by_support_threshold(sample, options.support_threshold);

    // write results of current last stage aggregator node to a part file
    std::string part_file_name(outdir + "epa." + std::to_string(local_rank) 
      + "." + std::to_string(chunk_num) + ".part");
    std::ofstream part_file(part_file_name);
    part_file << sample_to_jplace_string(sample, msa_stream);
    part_names.push_back(part_file_name);
    part_file.close();

#ifdef __MPI
    } // endif aggregate cleanup

    //==============================================================
    // EPA_MPI_STAGE_2_AGGREGATE === END
    //==============================================================
    timer.stop(); // stop timer of any stage

    if ( !(chunk_num % rebalance) ) // time to rebalance
    {
      lgr.dbg() << "Rebalancing..." << std::endl;
      int foreman = schedule[local_stage][0];
      // Step 1: aggregate the runtime statistics, first at the lowest rank per stage
      lgr.dbg() << "aggregate the runtime statistics..." << std::endl;
      epa_mpi_gather(timer, foreman, schedule[local_stage], local_rank);
      lgr.dbg() << "Runtime aggregate done!" << std::endl;
      
      // Step 2: calculate average time needed per chunk for the stage
      std::vector<double> perstage_avg(num_stages);

      int color = (local_rank == foreman) ? 1 : MPI_UNDEFINED;
      MPI_Comm foreman_comm;
      MPI_Comm_split(MPI_COMM_WORLD, color, local_rank, &foreman_comm);

      if (local_rank == foreman)
      {
        double avg = timer.average();
        // Step 3: make known to all other stage representatives (mpi_allgather)
        lgr.dbg() << "Foremen allgather..." << std::endl;
        MPI_Allgather(&avg, 1, MPI_DOUBLE, &perstage_avg[0], 1, MPI_DOUBLE, foreman_comm);
        lgr.dbg() << "Foremen allgather done!" << std::endl;
        MPI_Comm_free(&foreman_comm);
      }
      MPI_BARRIER(MPI_COMM_WORLD);
      // Step 4: stage representatives forward results to all stage members
      // epa_mpi_bcast(perstage_avg, foreman, schedule[local_stage], local_rank);
      lgr.dbg() << "Broadcasting..." << std::endl;
      MPI_Comm stage_comm;
      MPI_Comm_split(MPI_COMM_WORLD, local_stage, local_rank, &stage_comm);
      MPI_Bcast(&perstage_avg[0], num_stages, MPI_DOUBLE, 0, stage_comm);
      MPI_Comm_free(&stage_comm);
      lgr.dbg() << "Broadcasting done!" << std::endl;

      // Step 5: calculate schedule on every rank, deterministically!
      to_difficulty(perstage_avg);
      auto nps = solve(num_stages, world_size, perstage_avg);
      reassign(local_rank, nps, schedule, &local_stage);
      // Step 6: re-engage pipeline with new assignments
      // compute stages should try to keep their edge assignment! affinity!
      lgr.dbg() << "Rebalancing done!" << std::endl;
    }

#endif // __MPI
    sample.clear();
    msa_stream.clear();
    lgr.dbg() << "Chunk " << chunk_num << " done!" << std::endl;
    chunk_num++;
  }


  //==============================================================
  // POST COMPUTATION
  //==============================================================
  lgr.dbg() << "Starting Post-Comp" << std::endl;
  // finally, paste all part files together
#ifdef __MPI
  MPI_BARRIER(MPI_COMM_WORLD);

  if (local_rank != EPA_MPI_DEDICATED_WRITE_RANK)
  {
    if (local_stage == EPA_MPI_STAGE_LAST_AGGREGATE)
      epa_mpi_send(part_names, EPA_MPI_DEDICATED_WRITE_RANK, MPI_COMM_WORLD);
  }
  else
  {
    for (auto rank : schedule[EPA_MPI_STAGE_LAST_AGGREGATE])
    {
      if (rank == local_rank)
        continue;
      std::vector<std::string> remote_obj;
      epa_mpi_recieve(remote_obj, rank, MPI_COMM_WORLD);
      // merge(part_names, remote_obj);
      part_names.insert(part_names.end(), remote_obj.begin(), remote_obj.end());
    }
#endif
    // create output file
    std::ofstream outfile(outdir + "epa_result.jplace");
    lgr << "\nOutput file: " << outdir + "epa_result.jplace" << std::endl;
    outfile << init_jplace_string(get_numbered_newick_string(reference_tree.tree()));
    merge_into(outfile, part_names);
    outfile << finalize_jplace_string(invocation);
    outfile.close();
#ifdef __MPI
  }
#endif
}

// ================== LEGACY CODE ==========================================


Sample place(Tree& reference_tree, MSA& query_msa_)
{
  auto options = reference_tree.options();

  const auto num_branches = reference_tree.nums().branches;
  const auto num_queries = query_msa_.size();
  // get all edges
  std::vector<pll_utree_t *> branches(num_branches);
  auto num_traversed_branches = utree_query_branches(reference_tree.tree(), &branches[0]);
  assert(num_traversed_branches == num_branches);

  lgr << "\nPlacing "<< std::to_string(num_queries) << " reads on " <<
    std::to_string(num_branches) << " branches." << std::endl;

  // build all tiny trees with corresponding edges
  std::vector<Tiny_Tree> insertion_trees;
  for (size_t branch_id = 0; branch_id < num_branches; ++branch_id)
    insertion_trees.emplace_back(branches[branch_id], branch_id, reference_tree, !options.prescoring);
    /* clarification: last arg here is a flag specifying whether to optimize the branches.
      we don't want that if the mode is prescoring */

  // output class
  Sample sample(get_numbered_newick_string(reference_tree.tree()));
  for (size_t sequence_id = 0; sequence_id < num_queries; sequence_id++)
    sample.emplace_back(sequence_id, num_branches);


  // place all s on every edge
#ifdef __OMP
  #pragma omp parallel for schedule(dynamic)
#endif 
  for (size_t branch_id = 0; branch_id < num_branches; ++branch_id)
  {
    for (size_t sequence_id = 0; sequence_id < num_queries; ++sequence_id)
    {
      sample[sequence_id][branch_id] = insertion_trees[branch_id].place(query_msa_[sequence_id]);
      // printf("sequence %d branch %d:%f\n",sequence_id, branch_id, sample[sequence_id][branch_id].likelihood());
    }
  }
  // now that everything has been placed, we can compute the likelihood weight ratio
  compute_and_set_lwr(sample);

  /* prescoring was chosen: perform a second round, but only on candidate edges identified
    during the first run */
  if (options.prescoring)
  {
    lgr << "Entering second phase of placement. \n";
    if (options.prescoring_by_percentage)
      discard_bottom_x_percent(sample, (1.0 - options.prescoring_threshold));
    else
      discard_by_accumulated_threshold(sample, options.prescoring_threshold);

    // build a list of placements per edge that need to be recomputed
    std::vector<std::vector<std::tuple<Placement *, const unsigned int>>> recompute_list(num_branches);
    for (auto & pq : sample)
      for (auto & placement : pq)
        recompute_list[placement.branch_id()].push_back(std::make_tuple(&placement, pq.sequence_id()));

#ifdef __OMP
    #pragma omp parallel for schedule(dynamic)
#endif 
    for (size_t branch_id = 0; branch_id < num_branches; branch_id++)
    {
      Placement * placement;
      // Sequence * sequence;
      auto& branch = insertion_trees[branch_id];
      branch.opt_branches(true); // TODO only needs to be done once
      for (auto recomp_tuple : recompute_list[branch_id])
      {
        placement = std::get<0>(recomp_tuple);
        *placement = branch.place(query_msa_[std::get<1>(recomp_tuple)]);
      }
    }
    compute_and_set_lwr(sample);
  }

  // finally, trim the output
  if (options.acc_threshold)
    discard_by_accumulated_threshold(sample, options.support_threshold);
  else
    discard_by_support_threshold(sample, options.support_threshold);

  return sample;
}
