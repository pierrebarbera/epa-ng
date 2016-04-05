#include "place.hpp"

#include <stdexcept>
#include <fstream>
#include <chrono>

#include "file_io.hpp"
#include "jplace_util.hpp"
#include "stringify.hpp"
#include "set_manipulators.hpp"
#include "logging.hpp"
#include "Tiny_Tree.hpp"
#include "mpihead.hpp"

using namespace std;

Log lgr;

/**
  Schedules MPI processes to pipeline stages
*/
static int assign_color(int rank, int )
{
  int color = 0;
  if (rank != 0)
    color = EPA_MPI_STAGE_1_COMPUTE;
  else
    color = EPA_MPI_STAGE_1_AGGREGATE;

  return color;
}

void process(Tree& tree, MSA_Stream& msa_stream, const string& outdir,
              const Options& options, const string& invocation)
{

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

  unsigned int num_branches = 0;
  unsigned int chunk_size = 100;
  unsigned int num_sequences;

  // build all tiny trees with corresponding edges
  vector<Tiny_Tree> insertion_trees;
  for (auto node : branches)
    insertion_trees.emplace_back(node, partition, model, !options.prescoring);

  // create output file
  ofstream outfile(outdir + "epa_result.jplace");
  lgr << "\nOutput file: " << outdir + "epa_result.jplace" << endl;
  outfile << init_jplace(get_numbered_newick_string(tree));

  // output class
  Sample sample();

  // while data
  // TODO this could be a stream read, such that cat msa.fasta | epa .... would work
  while ((num_sequences = msa.read_next(chunk_size)) > 0)
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
        sample[cur_seq_id][local_branch_id] = insertion_trees[local_branch_id].place(msa[cur_seq_id]);
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

    if (options_.acc_threshold)
      discard_by_accumulated_threshold(sample, options_.support_threshold);
    else
      discard_by_support_threshold(sample, options_.support_threshold);

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
    msa.clear();
  }
  outfile << finalize_jplace(invocation);
  outfile.close();
}
