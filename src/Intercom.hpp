#pragma once

#include "schedule.hpp"
#include "mpihead.hpp"
#include "epa_mpi_util.hpp"
#include "Log.hpp"
#include "Timer.hpp"

#ifdef __MPI
#include <mpi.h>


/**
 * Class encapsulating the Communication and scheduling work used by the
 * distributed pipeline.
 */
class Intercom
{
public:
  Intercom(const size_t num_stages)
  {
    std::vector<double> initial_difficulty(num_stages, 1.0);
    MPI_Comm_rank(MPI_COMM_WORLD, &local_rank_);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size_);

    // get initial schedule
    auto init_nps = solve(num_stages, world_size_, initial_difficulty);
    assign(local_rank_, init_nps, schedule_, &local_stage_);

    lgr.dbg() << "Schedule: ";
    for (size_t i = 0; i < schedule_.size(); ++i) {
      lgr.dbg() << schedule_[i].size() << " ";
    }
    lgr.dbg() << std::endl;
  }

  Intercom()   = delete;
  ~Intercom()  = default;

  /**
   * Returns a set of MPI Ranks representing a stage in the pipeline.
   * 
   * TODO: may have to do a mapping from Typed_Stage::id to internal MPI stage id
   * in the future
   * 
   * @param  id index of a stage
   * @return    reference to a set of mpi ranks representing a stage
   */
  auto& schedule(const size_t id)
  {
    return schedule_[id];
  }

  /**
   * Returns a reference to the previously started async send requests.
   * 
   * TODO: eventually, to allow multiple independent, discontiguous MPI-stages
   *       per MPI rank, this should be indexed by the stage for which the
   *       previous requests are relevant
   *       
   * @return reference to structure holding previous async send requests information
   */
  auto& previous_requests()
  {
    return prev_requests_;
  }

  /**
   * Returns true when the supplied stage id should be active on the current MPI Rank.
   *
   * TODO: when mapping will be required, this functions impl. will change
   * 
   * @param  stage_id the ID of the stage
   * @return          wether that stage ought to be active on the current MPI Rank
   */
  bool stage_active(const size_t stage_id) const
  {
    return stage_id == static_cast<size_t>(local_stage_);
  }

  /**
   * Calculates, and propagates, a global pipeline schedule, based on local timing 
   * measurements. 
   * 
   * Synchronizes all MPI ranks!
   * 
   * @param timer local timing value used to calcuate the difficulty of a stage, and
   *              consequently the global schedule.
   */
  void rebalance(Timer& timer) 
  {
    MPI_BARRIER(MPI_COMM_WORLD);
    // if (local_rank == 0) 
    // {
    //   flight_time.stop();
    //   double aft = flight_time.average() / rebalance_delta;
    //   flight_file << chunk_num << ";" << aft << std::endl;
    // }

    lgr.dbg() << "Rebalancing..." << std::endl;
    const auto foreman = schedule_[local_stage_][0];
    const auto num_stages = schedule_.size(); 
    // Step 0: get per node average
    Timer per_node_avg({timer.average()});
    // Step 1: aggregate the runtime statistics, first at the lowest rank per stage
    lgr.dbg() << "aggregate the runtime statistics..." << std::endl;
    Timer dummy;
    epa_mpi_gather(per_node_avg, foreman, schedule_[local_stage_], local_rank_, dummy);
    lgr.dbg() << "Runtime aggregate done!" << std::endl;

    // Step 2: calculate total time needed per chunk for the stage, reflecting effort spent
    std::vector<double> perstage_total(num_stages);

    int color = (local_rank_ == foreman) ? 1 : MPI_UNDEFINED;
    MPI_Comm foreman_comm;
    MPI_Comm_split(MPI_COMM_WORLD, color, local_stage_, &foreman_comm);

    if (local_rank_ == foreman)
    {
      auto total_stagetime = per_node_avg.sum();
      // Step 3: make known to all other stage representatives (mpi_allgather)
      lgr.dbg() << "Foremen allgather..." << std::endl;
      MPI_Allgather(&total_stagetime, 1, MPI_DOUBLE, &perstage_total[0], 1, MPI_DOUBLE, foreman_comm);
      lgr.dbg() << "Foremen allgather done!" << std::endl;
      MPI_Comm_free(&foreman_comm);
    }
    // ensure all messages were received and previous requests are cleared
    epa_mpi_waitall(prev_requests_);
    
    MPI_BARRIER(MPI_COMM_WORLD);
    // Step 4: stage representatives forward results to all stage members
    // epa_mpi_bcast(perstage_total, foreman, schedule_[local_stage_], local_rank_);
    lgr.dbg() << "Broadcasting..." << std::endl;
    MPI_Comm stage_comm;
    MPI_Comm_split(MPI_COMM_WORLD, local_stage_, local_rank_, &stage_comm);

    // foreman must always be rank 0 in stage communicator:
    int split_key = local_rank_ == foreman ? -1 : local_rank_;
    MPI_Comm_split(MPI_COMM_WORLD, local_stage_, split_key, &stage_comm);
    MPI_Bcast(&perstage_total[0], num_stages, MPI_DOUBLE, 0, stage_comm);
    MPI_Comm_free(&stage_comm);
    lgr.dbg() << "Broadcasting done!" << std::endl;

    lgr.dbg() << "perstage total:";
    for (size_t i = 0; i < perstage_total.size(); ++i)
    {
      lgr.dbg() << " " << perstage_total[i];
    }
    lgr.dbg() << std::endl;

    // Step 5: calculate schedule_ on every rank, deterministically!
    to_difficulty(perstage_total);

    lgr.dbg() << "perstage difficulty:";
    for (size_t i = 0; i < perstage_total.size(); ++i)
    {
      lgr.dbg() << " " << perstage_total[i];
    }
    lgr.dbg() << std::endl;

    auto nps = solve(num_stages, world_size_, perstage_total);
    reassign(local_rank_, nps, schedule_, &local_stage_);
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
    // rebalance_delta *= 2;
    // rebalance += rebalance_delta;
    // reassign_happened = true;
    // if (local_rank_ == 0) 
    // {
    //   flight_time.clear();
    //   flight_time.start();
    // }

    prev_requests_.clear();
    timer.clear();
  }

  void barrier() const
  {
    MPI_BARRIER(MPI_COMM_WORLD);
  }

  int rank()
  {
    return local_rank_;
  }
 
private:
  int local_rank_   = -1;
  int world_size_   = -1;
  int local_stage_  = -1;

  schedule_type schedule_;
  previous_request_storage_t prev_requests_;

};

#else

class Intercom
{
public:
  Intercom(const size_t) {}
  Intercom() = default;
  ~Intercom()= default;


  // auto& schedule(const size_t) { }
  // auto& previous_requests() { }
  bool stage_active(const size_t) const { return true; }
  void rebalance(Timer&) { } 
  void barrier() const { }
  int rank() { return 0; }
  
};

#endif //__MPI
