#include "place.hpp"

#include <fstream>
#include <chrono>
#include <omp.h>

#include "file_io.hpp"
#include "jplace_util.hpp"
#include "stringify.hpp"
#include "set_manipulators.hpp"
#include "logging.hpp"
#include "Tiny_Tree.hpp"
#include "mpihead.hpp"
#include "pll_util.hpp"
#include "epa_pll_util.hpp"
#include "Timer.hpp"
#include "schedule.hpp"

#ifdef __MPI
#include "epa_mpi_util.hpp"
#endif

void process(Tree& reference_tree, MSA_Stream& msa_stream, const std::string& outdir,
              const Options& options, const std::string& invocation)
{
  auto model = reference_tree.model();

  int local_rank = 0;
#ifdef __MPI
  int world_size;
  MPI_Comm_rank(MPI_COMM_WORLD, &local_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);

  // shuffle available nodes to the different stages
  std::unordered_map<int, std::unordered_map<int, int>> global_rank;
  std::unordered_map<int, int> stage_size;
  int stage_1_compute_size = 0;
  int stage_1_aggregate_size = 0;
  int stage_2_compute_size = 0;
  int stage_2_aggregate_size = 0;
  int local_stage;
  const unsigned int num_stages = options.prescoring ? 4 : 2;

  unsigned int rebalance = 10;


  for (int i = 0; i < world_size; i++)
  {
    if (i == EPA_MPI_STAGE_1_AGGREGATE)
    {
      global_rank[EPA_MPI_STAGE_1_AGGREGATE][stage_1_aggregate_size++] = i;
      if (local_rank == i)
        local_stage = EPA_MPI_STAGE_1_AGGREGATE;
    }
    else if (i == EPA_MPI_STAGE_1_COMPUTE)
    {
      global_rank[EPA_MPI_STAGE_1_COMPUTE][stage_1_compute_size++] = i;
      if (local_rank == i)
        local_stage = EPA_MPI_STAGE_1_COMPUTE;
    }
    else if (i == EPA_MPI_STAGE_2_AGGREGATE)
    {
      global_rank[EPA_MPI_STAGE_2_AGGREGATE][stage_2_aggregate_size++] = i;
      if (local_rank == i)
        local_stage = EPA_MPI_STAGE_2_AGGREGATE;
    }
    else if (i == EPA_MPI_STAGE_2_COMPUTE)
    {
      global_rank[EPA_MPI_STAGE_2_COMPUTE][stage_2_compute_size++] = i;
      if (local_rank == i)
        local_stage = EPA_MPI_STAGE_2_COMPUTE;
    }
  }
  stage_size[EPA_MPI_STAGE_1_COMPUTE] = stage_1_compute_size;
  stage_size[EPA_MPI_STAGE_1_AGGREGATE] = stage_1_aggregate_size;
  stage_size[EPA_MPI_STAGE_2_COMPUTE] = stage_2_compute_size;
  stage_size[EPA_MPI_STAGE_2_AGGREGATE] = stage_2_aggregate_size;

  Timer timer;

  const auto EPA_MPI_STAGE_LAST_AGGREGATE 
    = options.prescoring ? EPA_MPI_STAGE_2_AGGREGATE : EPA_MPI_STAGE_1_AGGREGATE;
  const auto EPA_MPI_DEDICATED_WRITE_RANK = 0; // TODO should really be one of the last rank aggr

  if (local_rank == 0)
  {
#endif // __MPI
  lgr << "EPA - Evolutionary Placement Algorithm\n";
  lgr << "\nInvocation: \n" << invocation << "\n";
#ifdef __MPI
  } else {
    std::cout << "Rank " << local_rank << " checking in." << std::endl;
  }
#endif // __MPI

  // prepare all structures
  // TODO build edge set / OOC tree

  unsigned int chunk_size = 100;
  unsigned int num_sequences;

  const auto num_branches = reference_tree.nums().branches;

  // get all edges
  std::vector<pll_utree_t *> branches(num_branches);
  auto num_traversed_branches = utree_query_branches(reference_tree.tree(), &branches[0]);
  assert(num_traversed_branches == num_branches);

  // build all tiny trees with corresponding edges
  std::vector<Tiny_Tree> insertion_trees;
  for (unsigned int branch_id = 0; branch_id < num_branches; branch_id++)
  {
    // TODO check if current mpi node is supposed to get this branch id
    insertion_trees.emplace_back(branches[branch_id], branch_id, reference_tree, !options.prescoring);
  }

  // output class
  Sample sample;

  auto chunk_num = 0;
  std::vector<std::string> part_names;


  // while data
  // TODO this could be a stream read, such that cat msa.fasta | epa ... would work
  while ((num_sequences = msa_stream.read_next(chunk_size)) > 0)
  {
#ifdef __MPI
    timer.start(); // start timer of any stage
    //==============================================================
    // EPA_MPI_STAGE_1_COMPUTE === BEGIN
    //==============================================================
    if (local_stage == EPA_MPI_STAGE_1_COMPUTE)
    {
#endif // __MPI
    // prepare placement structure
    for (unsigned int cur_seq_id = 0; cur_seq_id < num_sequences; cur_seq_id++)
      sample.emplace_back(cur_seq_id, num_branches);

    // place sequences
    #pragma omp parallel for schedule(dynamic)
    for (unsigned int local_branch_id = 0; local_branch_id < num_branches; ++local_branch_id)
    {
      for (unsigned int cur_seq_id = 0; cur_seq_id < num_sequences; cur_seq_id++)
      {
        sample[cur_seq_id][local_branch_id] = insertion_trees[local_branch_id].place(msa_stream[cur_seq_id]);
      }
    }

#ifdef __MPI
    // MPI: split the result and send the part to correct aggregate node
    std::vector<Sample> parts;
    split(sample, parts, global_rank[EPA_MPI_STAGE_1_AGGREGATE].size());
    for (unsigned int i = 0; i < parts.size(); ++i)
    {
      epa_mpi_send(parts[i], global_rank[EPA_MPI_STAGE_1_AGGREGATE][i], MPI_COMM_WORLD);
    }

    // MPI Barrier first compute stage TODO for now implicit barrier by synchronous comm

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
    for (auto rank_pair : global_rank[EPA_MPI_STAGE_1_COMPUTE])
    {
      int rank = rank_pair.second;
      Sample remote_sample;
      epa_mpi_recieve(remote_sample, rank, MPI_COMM_WORLD);
      merge(sample, remote_sample);
    }
#endif // __MPI
    // build lwrs
    compute_and_set_lwr(sample);

    // if this was a prescring run, select the candidate edges
    if(options.prescoring)
    {
      if (options.prescoring_by_percentage)
        discard_bottom_x_percent(sample, (1.0 - options.prescoring_threshold));
      else
        discard_by_accumulated_threshold(sample, options.prescoring_threshold);
    }
#ifdef __MPI
    // if prescoring was selected, we need to send the intermediate results off to thorough placement
    if(options.prescoring)
    {
      // TODO split!
      epa_mpi_send(sample, global_rank[EPA_MPI_STAGE_2_COMPUTE][0], MPI_COMM_WORLD);
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
    // (MPI: recieve results, merge them)
    for (auto rank_pair : global_rank[EPA_MPI_STAGE_1_AGGREGATE])
    {
      int rank = rank_pair.second;
      Sample remote_sample; 
      epa_mpi_recieve(remote_sample, rank, MPI_COMM_WORLD);
      merge(sample, remote_sample);
    }
#endif // __MPI
    if (options.prescoring)
    {
      // build a list of placements per edge that need to be recomputed
      std::vector<std::vector<std::tuple<Placement *, const unsigned int>>> recompute_list(num_branches);
      for (auto & pq : sample)
        for (auto & placement : pq)
          recompute_list[placement.branch_id()].push_back(std::make_tuple(&placement, pq.sequence_id()));

      #pragma omp parallel for schedule(dynamic)
      for (unsigned int branch_id = 0; branch_id < num_branches; branch_id++)
      {
        Placement * placement;
        // Sequence * sequence;
        auto& branch = insertion_trees[branch_id];
        branch.opt_branches(true); // TODO only needs to be done once
        for (auto recomp_tuple : recompute_list[branch_id])
        {
          placement = std::get<0>(recomp_tuple);
          *placement = branch.place(msa_stream[std::get<1>(recomp_tuple)]);
        }
      }
    }
#ifdef __MPI
    if(options.prescoring)
    {
      std::vector<Sample> parts;
      split(sample, parts, global_rank[EPA_MPI_STAGE_2_AGGREGATE].size());
      for (unsigned int i = 0; i < parts.size(); ++i)
      {
        epa_mpi_send(parts[i], global_rank[EPA_MPI_STAGE_2_AGGREGATE][i], MPI_COMM_WORLD);
      }
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
    // (MPI: recieve results, merge them)
    for (auto rank_pair : global_rank[EPA_MPI_STAGE_2_COMPUTE])
    {
      int rank = rank_pair.second;
      Sample remote_sample;
      epa_mpi_recieve(remote_sample, rank, MPI_COMM_WORLD);
      merge(sample, remote_sample);
    }
    }
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
    timer.stop(); // stop timer of any stage

    if ( !(chunk_num % rebalance) ) // time to rebalance
    {
      int foreman = global_rank[local_stage][0];
      // Step 1: aggregate the runtime statistics, first at the lowest rank per stage
      epa_mpi_gather(timer, foreman, global_rank[local_stage], local_rank);
      
      // Step 2: calculate average time needed per chunk for the stage
      int color = MPI_UNDEFINED;
      std::vector<double> perstage_avg(num_stages);

      if (local_rank == foreman)
      {
        color = 1;
        double avg = timer.average();
        // Step 3: make known to all other stage representatives (mpi_allgather)
        MPI_Comm foreman_comm;
        MPI_Comm_split(MPI_COMM_WORLD, color, local_rank, &foreman_comm);
        MPI_Allgather(&avg, 1, MPI_DOUBLE, &perstage_avg[0], 1, MPI_DOUBLE, foreman_comm);
        MPI_Comm_free(&foreman_comm);
      }
      // Step 4: stage representatives forward results to all stage members
      // epa_mpi_bcast(perstage_avg, foreman, global_rank[local_stage], local_rank);
      MPI_Comm stage_comm;
      MPI_Comm_split(MPI_COMM_WORLD, local_stage, local_rank, &stage_comm);
      MPI_Bcast(&perstage_avg[0], num_stages, MPI_DOUBLE, 0, stage_comm);
      MPI_Comm_free(&stage_comm);

      // Step 5: calculate schedule on every rank, deterministically!
      to_difficulty(perstage_avg);
      auto sched = solve(num_stages, world_size, perstage_avg);
      // Step 6: re-engage pipeline with new assignments
      // compute stages should try to keep their edge assignment! affinity!
    }

#endif // __MPI
    sample.clear();
    msa_stream.clear();
    chunk_num++;
  }


  //==============================================================
  // POST COMPUTATION
  //==============================================================
  // finally, paste all part files together
#ifdef __MPI
  MPI_BARRIER(MPI_COMM_WORLD);

  if (local_rank != EPA_MPI_DEDICATED_WRITE_RANK)
  {
    if (local_stage == EPA_MPI_STAGE_LAST_AGGREGATE)
    {
      epa_mpi_send(part_names, EPA_MPI_DEDICATED_WRITE_RANK, MPI_COMM_WORLD);
    }
  }
  else
  {
  for (auto rank_pair : global_rank[EPA_MPI_STAGE_LAST_AGGREGATE])
  {
    int rank = rank_pair.second;
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
  for (unsigned int branch_id = 0; branch_id < num_branches; ++branch_id)
    insertion_trees.emplace_back(branches[branch_id], branch_id, reference_tree, !options.prescoring);
    /* clarification: last arg here is a flag specifying whether to optimize the branches.
      we don't want that if the mode is prescoring */

  // output class
  Sample sample(get_numbered_newick_string(reference_tree.tree()));
  for (unsigned int sequence_id = 0; sequence_id < num_queries; sequence_id++)
    sample.emplace_back(sequence_id, num_branches);


  // place all s on every edge
  #pragma omp parallel for schedule(dynamic)
  for (unsigned int branch_id = 0; branch_id < num_branches; ++branch_id)
  {
    for (unsigned int sequence_id = 0; sequence_id < num_queries; ++sequence_id)
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

    #pragma omp parallel for schedule(dynamic)
    for (unsigned int branch_id = 0; branch_id < num_branches; branch_id++)
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
