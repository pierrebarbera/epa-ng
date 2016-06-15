#include "place.hpp"

#include <stdexcept>
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

#ifdef __MPI
#include "epa_mpi_util.hpp"
#endif

using namespace std;

Log lgr;

void process(Tree& reference_tree, MSA_Stream& msa_stream, const string& outdir,
              const Options& options, const string& invocation)
{
  auto model = reference_tree.model();

  lgr = Log(outdir + "epa_info.log");

#ifdef __MPI
  int world_rank, world_size;
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);

  // shuffle available nodes to the different stages
  unordered_map<int, unordered_map<int, int>> global_rank;
  unordered_map<int, int> stage_size;
  int stage_1_compute_size = 0;
  int stage_1_aggregate_size = 0;
  int stage_2_compute_size = 0;
  int stage_2_aggregate_size = 0;
  int local_stage;


  for (int i = 0; i < world_size; i++)
  {
    if (i == 0)
    {
      global_rank[EPA_MPI_STAGE_1_AGGREGATE][stage_1_aggregate_size++] = i;
      if (world_rank == i)
        local_stage = EPA_MPI_STAGE_1_AGGREGATE;
    }
    else if (i == 1)
    {
      global_rank[EPA_MPI_STAGE_1_COMPUTE][stage_1_compute_size++] = i;
      if (world_rank == i)
        local_stage = EPA_MPI_STAGE_1_COMPUTE;
    }
    else if (i == 2)
    {
      global_rank[EPA_MPI_STAGE_2_AGGREGATE][stage_2_aggregate_size++] = i;
      if (world_rank == i)
        local_stage = EPA_MPI_STAGE_2_AGGREGATE;
    }
    else if (i == 3)
    {
      global_rank[EPA_MPI_STAGE_2_COMPUTE][stage_2_compute_size++] = i;
      if (world_rank == i)
        local_stage = EPA_MPI_STAGE_2_COMPUTE;
    }
  }
  stage_size[EPA_MPI_STAGE_1_COMPUTE] = stage_1_compute_size;
  stage_size[EPA_MPI_STAGE_1_AGGREGATE] = stage_1_aggregate_size;
  stage_size[EPA_MPI_STAGE_2_COMPUTE] = stage_2_compute_size;
  stage_size[EPA_MPI_STAGE_2_AGGREGATE] = stage_2_aggregate_size;

  Timer timer;

  if (world_rank == 0)
  {
#endif // __MPI
  lgr << "EPA - Evolutionary Placement Algorithm\n";
  lgr << "\nInvocation: \n" << invocation << "\n";
#ifdef __MPI
  } else {
    lgr << "Rank " << world_rank << " checking in.\n";
  }
#endif // __MPI

  // prepare all structures
  // TODO build edge set / OOC tree

  unsigned int chunk_size = 100;
  unsigned int num_sequences;

  const auto num_branches = reference_tree.nums().branches;

  // get all edges
  vector<pll_utree_t *> branches(num_branches);
  auto num_traversed_branches = utree_query_branches(reference_tree.tree(), &branches[0]);
  assert(num_traversed_branches == num_branches);

  // build all tiny trees with corresponding edges
  vector<Tiny_Tree> insertion_trees;
  for (unsigned int branch_id = 0; branch_id < num_branches; branch_id++)
  {
    // TODO check if current mpi node is supposed to get this branch id
    insertion_trees.emplace_back(branches[branch_id], branch_id, reference_tree, !options.prescoring);
  }

  // create output file
  ofstream outfile(outdir + "epa_result.jplace");
  lgr << "\nOutput file: " << outdir + "epa_result.jplace" << endl;
  outfile << init_jplace_string(get_numbered_newick_string(reference_tree.tree()));

  // output class
  Sample sample;

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
    // TODO split!
    epa_mpi_send(sample, global_rank[EPA_MPI_STAGE_1_AGGREGATE][0], MPI_COMM_WORLD);

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
    for (auto rank : global_rank[EPA_MPI_STAGE_1_COMPUTE])
    {
      (void)rank;// surpress warning
      Sample remote_sample;
      epa_mpi_recieve(remote_sample, MPI_ANY_SOURCE, MPI_COMM_WORLD);
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
    if (local_stage == EPA_MPI_STAGE_2_COMPUTE)
    {
    // (MPI: recieve results, merge them)
    for (auto rank : global_rank[EPA_MPI_STAGE_1_AGGREGATE])
    {
      (void)rank;// surpress warning
      Sample remote_sample;
      epa_mpi_recieve(remote_sample, MPI_ANY_SOURCE, MPI_COMM_WORLD);
      merge(sample, remote_sample);
    }
#endif // __MPI
    if (options.prescoring)
    {
      // build a list of placements per edge that need to be recomputed
      vector<vector<tuple<Placement *, const unsigned int>>> recompute_list(num_branches);
      for (auto & pq : sample)
        for (auto & placement : pq)
          recompute_list[placement.branch_id()].push_back(make_tuple(&placement, pq.sequence_id()));

      #pragma omp parallel for schedule(dynamic)
      for (unsigned int branch_id = 0; branch_id < num_branches; branch_id++)
      {
        Placement * placement;
        // Sequence * sequence;
        auto& branch = insertion_trees[branch_id];
        branch.opt_branches(true); // TODO only needs to be done once
        for (auto recomp_tuple : recompute_list[branch_id])
        {
          placement = get<0>(recomp_tuple);
          *placement = branch.place(msa_stream[get<1>(recomp_tuple)]);
        }
      }
    }
#ifdef __MPI
    if(options.prescoring)
    {
      // TODO split!
      epa_mpi_send(sample, global_rank[EPA_MPI_STAGE_2_AGGREGATE][0], MPI_COMM_WORLD);
    }
    } // endif (local_stage == EPA_MPI_STAGE_2_COMPUTE)
    //==============================================================
    // EPA_MPI_STAGE_2_COMPUTE === END
    //==============================================================

    //==============================================================
    // EPA_MPI_STAGE_2_AGGREGATE === BEGIN
    //==============================================================
    if ((local_stage == EPA_MPI_STAGE_1_AGGREGATE and not options.prescoring) or
        (local_stage == EPA_MPI_STAGE_2_AGGREGATE and options.prescoring))
    {
    // only if this is the 4th stage do we need to get from mpi
    if (local_stage == EPA_MPI_STAGE_2_AGGREGATE)
    {
    // (MPI: recieve results, merge them)
    for (auto rank : global_rank[EPA_MPI_STAGE_2_COMPUTE])
    {
      (void)rank;// surpress warning
      Sample remote_sample;
      epa_mpi_recieve(remote_sample, MPI_ANY_SOURCE, MPI_COMM_WORLD);
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

    // finally, write to file
    outfile << sample_to_jplace_string(sample, msa_stream);
#ifdef __MPI
    } // endif aggregate cleanup
    timer.stop(); // stop timer of any stage
#endif // __MPI
    sample.clear();
    msa_stream.clear();
  }
  outfile << finalize_jplace_string(invocation);
  outfile.close();
}

// ================== LEGACY CODE ==========================================


Sample place(Tree& reference_tree, MSA& query_msa_)
{
  auto options = reference_tree.options();

  const auto num_branches = reference_tree.nums().branches;
  const auto num_queries = query_msa_.size();
  // get all edges
  vector<pll_utree_t *> branches(num_branches);
  auto num_traversed_branches = utree_query_branches(reference_tree.tree(), &branches[0]);
  assert(num_traversed_branches == num_branches);

  lgr << "\nPlacing "<< to_string(num_queries) << " reads on " <<
    to_string(num_branches) << " branches." << endl;

  // build all tiny trees with corresponding edges
  vector<Tiny_Tree> insertion_trees;
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
    vector<vector<tuple<Placement *, const unsigned int>>> recompute_list(num_branches);
    for (auto & pq : sample)
      for (auto & placement : pq)
        recompute_list[placement.branch_id()].push_back(make_tuple(&placement, pq.sequence_id()));

    #pragma omp parallel for schedule(dynamic)
    for (unsigned int branch_id = 0; branch_id < num_branches; branch_id++)
    {
      Placement * placement;
      // Sequence * sequence;
      auto& branch = insertion_trees[branch_id];
      branch.opt_branches(true); // TODO only needs to be done once
      for (auto recomp_tuple : recompute_list[branch_id])
      {
        placement = get<0>(recomp_tuple);
        *placement = branch.place(query_msa_[get<1>(recomp_tuple)]);

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
