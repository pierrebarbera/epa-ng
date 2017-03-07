#pragma once

#ifdef __MPI

#include "mpihead.hpp"
#include "Sample.hpp"
#include "Timer.hpp"

#include <sstream>
#include <memory>
#include <unordered_map>
#include <cereal/archives/binary.hpp>
#include <cereal/types/vector.hpp>
#include <cereal/types/string.hpp>

// types to keep track of previous async sends
typedef struct
{
  MPI_Request  req;
  char*         buf = nullptr;
} request_tuple_t;

typedef typename std::unordered_map<int, request_tuple_t>  previous_request_storage_t;

static void err_check(int errval)
{
  if (errval != MPI_SUCCESS)
  {
    std::stringstream msg;
    msg << "Failed MPI call: ";
    switch (errval) {
      case MPI_ERR_COMM:
        msg << "Invalid communicator.";
        break;
      case MPI_ERR_TYPE:
        msg << "Invalid datatype argument.";
        break;
      case MPI_ERR_COUNT:
        msg << "Invalid count argument.";
        break;
      case MPI_ERR_TAG:
        msg << "Invalid tag argument.";
        break;
      case MPI_ERR_RANK:
        msg << "Invalid source or destination rank.";
        break;
      default:
        msg << "unknown";
    }
    msg << std::endl;
    throw std::runtime_error{msg.str()};
  }
}

void epa_mpi_waitall(previous_request_storage_t reqs)
{
  for (auto& pair : reqs)
  {
    auto& r = pair.second;
    if (r.req)
    {
      MPI_Status status;
      err_check(MPI_Wait(&r.req, &status));
      delete[] r.buf;
      r.buf = nullptr;
    }
  }
}

template <typename T>
void epa_mpi_send(T& obj, int dest_rank, MPI_Comm comm)
{
  // serialize the obj
  std::stringstream ss;
  cereal::BinaryOutputArchive out_archive(ss);
  out_archive(obj);

  // send obj to specified node
  std::string data = ss.str();
  auto buffer = new char[data.size()];
  memcpy(buffer, data.c_str(), data.size() * sizeof(char));
  err_check(MPI_Send(buffer, data.size(), MPI_CHAR, dest_rank, 0, comm));
  delete[] buffer;
}

template <typename T>
void epa_mpi_isend(T& obj, int dest_rank, MPI_Comm comm, request_tuple_t& prev_req, Timer& timer)
{
  // wait for completion of previous send
  if (prev_req.req)
  {
    MPI_Status status;
    timer.pause();
    err_check(MPI_Wait(&prev_req.req, &status));
    timer.resume();
    delete[] prev_req.buf;
    // free previous request?
  }

  // serialize the obj
  std::stringstream ss;
  cereal::BinaryOutputArchive out_archive(ss);
  out_archive(obj);

  // send obj to specified node
  std::string data = ss.str();
  auto buffer = new char[data.size()];
  memcpy(buffer, data.c_str(), data.size() * sizeof(char));
  err_check(MPI_Issend(buffer, data.size(), MPI_CHAR, dest_rank, 0, comm, &prev_req.req));
  
  prev_req.buf = buffer;
}

template <typename T>
void epa_mpi_recieve(T& obj, int src_rank, MPI_Comm comm, Timer& timer)
{
  // probe to find out the message size
  MPI_Status status;
  int size;
  timer.pause();
  err_check(MPI_Probe(src_rank, 0, comm, &status));
  timer.resume();
  MPI_Get_count(&status, MPI_CHAR, &size);

  // prepare buffer
  auto buffer = new char[size];

  //  get the actual payload
  err_check(MPI_Recv(buffer, size, MPI_CHAR, src_rank, 0, comm, &status));

  // deserialization
  std::stringstream ss;
  ss.write(buffer, size);
  cereal::BinaryInputArchive in_archive(ss);

  // build the object
  in_archive(obj);

  delete[] buffer;
}

template <typename T>
static inline void isend_all(std::vector<T>& parts, std::vector<int>& dest_ranks, MPI_Comm comm, previous_request_storage_t& prev_reqs, Timer& timer)
{
  for (size_t i = 0; i < parts.size(); ++i)
  {
    auto dest = dest_ranks[i];
    epa_mpi_isend(parts[i], dest, comm, prev_reqs[dest], timer);
  }
}

template <typename T>
void epa_mpi_split_send(T& obj, std::vector<int>& dest_ranks, MPI_Comm comm, previous_request_storage_t& prev_reqs, Timer& timer)
{
  std::vector<T> parts;
  split(obj, parts, dest_ranks.size());
  isend_all(parts, dest_ranks, comm, prev_reqs, timer);
}

enum class recieve_status {WAITING, READY, DONE};

struct status_type
{
  int rank;
  recieve_status status;
};

// template <typename T>
// void epa_mpi_recieve_merge(T& obj, std::vector<int>& src_ranks, MPI_Comm comm, Timer& timer)
// {
//   const auto num_src_ranks = src_ranks.size();
//   size_t num_recieved = 0;

//   std::vector<status_type> status_list(num_src_ranks);
//   size_t i = 0;
//   for (auto& rank : src_ranks)
//   {
//     status_list[i].rank   = rank;
//     status_list[i].status = recieve_status::WAITING;
//     ++i;
//   }
  
//     to avoid waits as much as possible: build a list that tells us for which src ranks
//     there are currently messages waiting to be recieved, then recieve those.
//     Repeat this procedure of nonblocking check and blocking recieve until all source ranks
//     have been visited ONCE (as we still want per chunk synchronicity)
   
//   do
//   {
//     // nonblocking probe to check which are ready
//     for (auto& s : status_list)
//     {
//       MPI_Status stat;
//       int message_available = false;
//       err_check(MPI_Iprobe(s.rank, 0, comm, &message_available, &stat));
//       if (message_available) s.status = recieve_status::READY;
//     }

//     // blocking recieve of all that are ready
//     for (auto& s : status_list)
//     {
//       if (s.status == recieve_status::READY)
//       {
//         T remote_obj;
//         epa_mpi_recieve(remote_obj, s.rank, comm, timer);
//         merge(obj, remote_obj);
//         ++num_recieved;
//         s.status = recieve_status::DONE;
//       }
//     }
//   } while (num_recieved < num_src_ranks);
// }

template <typename T>
void epa_mpi_recieve_merge(T& obj, std::vector<int>& src_ranks, MPI_Comm comm, Timer& timer)
{
  for (auto rank : src_ranks)
  {
    T remote_obj;
    epa_mpi_recieve(remote_obj, rank, comm, timer);
    merge(obj, remote_obj);
  }
}
template <typename T>
void epa_mpi_gather(T& obj, int dest_rank, std::vector<int>& src_ranks, int local_rank, Timer& timer)
{
  if (dest_rank != local_rank)
  {
    epa_mpi_send(obj, dest_rank, MPI_COMM_WORLD);
  }
  else
  {
    for (auto src_rank : src_ranks)
    {
      if (local_rank == src_rank)
        continue;
      T remote_obj;
      epa_mpi_recieve(remote_obj, src_rank, MPI_COMM_WORLD, timer);
      merge(obj, remote_obj);
    }
  }
}

template <typename T>
void epa_mpi_bcast(T& obj, int src_rank, std::vector<int>& dest_ranks, int local_rank, Timer& timer)
{
  if (src_rank == local_rank)
  {
    for (auto dest_rank : dest_ranks)
    {
      if (local_rank == dest_rank)
        continue;
      epa_mpi_send(obj, dest_rank, MPI_COMM_WORLD);
    }
  }
  else
  {
    epa_mpi_recieve(obj, src_rank, MPI_COMM_WORLD, timer);
  }
}

#endif

