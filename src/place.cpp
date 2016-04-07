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

using namespace std;

Log lgr;

void process(Tree& epa_tree, MSA_Stream& msa_stream, const string& outdir,
              const Options& options, const string& invocation)
{
  auto model = epa_tree.model();
  auto partition = epa_tree.partition();
  auto tree = epa_tree.tree();

  lgr = Log(outdir + "epa_info.log");

  lgr << "EPA - Evolutionary Placement Algorithm" << endl;
  lgr << "\nInvocation: \n" << invocation << endl;


  #ifdef __MPI
  int stage_rank, stage_size, world_rank, world_size;
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);

  // shuffle available nodes to the different stages
  // map[stage_id][stage_rank] = global_rank
  unordered_map<int, unordered_map<int, int>> global_rank;
  // map[stage_id] = stage_size
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
    else
    {
      global_rank[EPA_MPI_STAGE_1_COMPUTE][stage_1_compute_size++] = i;
      if (world_rank == i)
        local_stage = EPA_MPI_STAGE_1_COMPUTE;
    }

  }
  stage_size[EPA_MPI_STAGE_1_COMPUTE] = stage_1_compute_size;
  stage_size[EPA_MPI_STAGE_1_AGGREGATE] = stage_1_aggregate_size;
  stage_size[EPA_MPI_STAGE_2_COMPUTE] = stage_2_compute_size;
  stage_size[EPA_MPI_STAGE_2_AGGREGATE] = stage_2_aggregate_size;

  #endif // __MPI

  // prepare all structures
  // TODO build edge set / OOC tree

  unsigned int chunk_size = 100;
  unsigned int num_sequences;

  auto nums = epa_tree.nums();
  const auto num_branches = nums.branches;

  // get all edges
  vector<pll_utree_t *> branches(num_branches);
  auto num_traversed_branches = utree_query_branches(tree, &branches[0]);
  assert(num_traversed_branches == num_branches);

  // build all tiny trees with corresponding edges
  vector<Tiny_Tree> insertion_trees;
  for (unsigned int branch_id = 0; branch_id < num_branches; branch_id++)
  {
    // TODO check if current mpi node is supposed to get this branch id
    auto node = branches[branch_id];
    insertion_trees.emplace_back(node, branch_id, partition, model, !options.prescoring);
  }

  // create output file
  ofstream outfile(outdir + "epa_result.jplace");
  lgr << "\nOutput file: " << outdir + "epa_result.jplace" << endl;
  outfile << init_jplace_string(get_numbered_newick_string(tree));

  // output class
  Sample sample;

  // while data
  // TODO this could be a stream read, such that cat msa.fasta | epa .... would work
  while ((num_sequences = msa_stream.read_next(chunk_size)) > 0)
  {
    #ifdef __MPI
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
    for (unsigned int local_branch_id = 0; local_branch_id < num_branches; ++local_branch_id)
    {
      for (unsigned int cur_seq_id = 0; cur_seq_id < num_sequences; cur_seq_id++)
      {
        sample[cur_seq_id][local_branch_id] = insertion_trees[local_branch_id].place(msa_stream[cur_seq_id]);
      }
    }

    #ifdef __MPI
    // MPI: split the result and send the part to correct aggregate node
    epa_mpi_send(sample, global_rank[EPA_MPI_STAGE_1_AGGREGATE][], MPI_COMM_WORLD);


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
      Sample remote_sample;
      epa_mpi_recieve(remote_sample, MPI_ANY_SOURCE, MPI_COMM_WORLD);
      merge(sample, remote_sample);
    }
    #endif // __MPI
    // build lwrs, discard
    compute_and_set_lwr(sample);

    if (options.acc_threshold)
      discard_by_accumulated_threshold(sample, options.support_threshold);
    else
      discard_by_support_threshold(sample, options.support_threshold);

    #ifdef __MPI
    } // endif (local_stage == EPA_MPI_STAGE_1_AGGREGATE)
    //==============================================================
    // EPA_MPI_STAGE_1_AGGREGATE === END
    //==============================================================
    #endif // __MPI


    // produce map of what sequence needs to recompute what
    // MPI: scatter all such partial maps to the second compute stage

    // take recompute list (MPI: recieve times number of aggregators)
    // MPI: compute what current rank needs to recompute
    // recompute based on list
    // produce placements (MPI: send to correct 2nd phase aggregator)

    // take re-placed placements (MPI: recieve)
    // rebuild their lwrs, do second stage discard
    // save to jplace (MPI: get file lock)

    sample.clear();
    msa_stream.clear();
  }
  outfile << finalize_jplace_string(invocation);
  outfile.close();
}

// ================== LEGACY CODE ==========================================


Sample place(Tree& epa_tree, MSA& query_msa_)
{
  auto options_ = epa_tree.options();
  auto model_ = epa_tree.model();
  auto tree_ = epa_tree.tree();
  auto partition_ = epa_tree.partition();
  auto nums_ = epa_tree.nums();

  const auto num_branches = nums_.branches;
  const auto num_queries = query_msa_.size();
  // get all edges
  vector<pll_utree_t *> branches(num_branches);
  auto num_traversed_branches = utree_query_branches(tree_, &branches[0]);
  assert(num_traversed_branches == num_branches);

  lgr << "\nPlacing "<< to_string(num_queries) << " reads on " <<
    to_string(num_branches) << " branches." << endl;

  // build all tiny trees with corresponding edges
  vector<Tiny_Tree> insertion_trees;
  for (unsigned int branch_id = 0; branch_id < num_branches; ++branch_id)
    insertion_trees.emplace_back(branches[branch_id], branch_id, partition_, model_, !options_.prescoring);
    /* clarification: last arg here is a flag specifying whether to optimize the branches.
      we don't want that if the mode is prescoring */

  // output class
  Sample sample(get_numbered_newick_string(tree_));
  for (unsigned int sequence_id = 0; sequence_id < num_queries; sequence_id++)
    sample.emplace_back(sequence_id, num_branches);


  // place all s on every edge
  #pragma omp parallel for schedule(dynamic)
  for (unsigned int branch_id = 0; branch_id < num_branches; ++branch_id)
  {
    for (unsigned int sequence_id = 0; sequence_id < num_queries; ++sequence_id)
    {
      sample[sequence_id][branch_id] = insertion_trees[branch_id].place(query_msa_[sequence_id]);
    }
  }
  // now that everything has been placed, we can compute the likelihood weight ratio
  compute_and_set_lwr(sample);

  /* prescoring was chosen: perform a second round, but only on candidate edges identified
    during the first run */
  if (options_.prescoring)
  {
    lgr << "Entering second phase of placement. \n";
    if (options_.prescoring_by_percentage)
      discard_bottom_x_percent(sample, (1.0 - options_.prescoring_threshold));
    else
      discard_by_accumulated_threshold(sample, options_.prescoring_threshold);

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
  if (options_.acc_threshold)
    discard_by_accumulated_threshold(sample, options_.support_threshold);
  else
    discard_by_support_threshold(sample, options_.support_threshold);

  return sample;
}
