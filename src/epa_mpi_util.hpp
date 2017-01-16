#pragma once

#ifdef __MPI

#include "mpihead.hpp"
#include "Sample.hpp"

#include <sstream>
#include <memory>
#include <unordered_map>
#include <cereal/archives/binary.hpp>
#include <cereal/types/vector.hpp>
#include <cereal/types/string.hpp>

// types to keep track of previous async sends
typedef struct
{
  MPI_Request*  req = nullptr;
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
      err_check(MPI_Wait(r.req, &status));
      delete[] r.buf;
      r.buf = nullptr;
      r.req = nullptr;
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
void epa_mpi_isend(T& obj, int dest_rank, MPI_Comm comm, request_tuple_t& prev_req)
{
  // wait for completion of previous send
  if (prev_req.req)
  {
    MPI_Status status;
    err_check(MPI_Wait(prev_req.req, &status));
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
  err_check(MPI_Isend(buffer, data.size(), MPI_CHAR, dest_rank, 0, comm, prev_req.req));
  
  prev_req.buf = buffer;
}

template <typename T>
void epa_mpi_recieve(T& obj, int src_rank, MPI_Comm comm)
{
  // probe to find out the message size
  MPI_Status status;
  int size;
  err_check(MPI_Probe(src_rank, 0, comm, &status));
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
static inline void isend_all(std::vector<T>& parts, std::vector<int>& dest_ranks, MPI_Comm comm, previous_request_storage_t& prev_reqs)
{
  for (size_t i = 0; i < parts.size(); ++i)
  {
    auto dest = dest_ranks[i];
    epa_mpi_isend(parts[i], dest, comm, prev_reqs[dest]);
  }
}

template <typename T>
void epa_mpi_split_send(T& obj, std::vector<int>& dest_ranks, MPI_Comm comm, previous_request_storage_t& prev_reqs)
{
  std::vector<T> parts;
  split(obj, parts, dest_ranks.size());
  isend_all(parts, dest_ranks, comm, prev_reqs);
}

void epa_mpi_split_send(Sample& obj, const unsigned int num_seq, std::vector<int>& dest_ranks, MPI_Comm comm, previous_request_storage_t& prev_reqs)
{
  std::vector<Sample> parts;
  split(obj, parts, dest_ranks.size(), num_seq);
  isend_all(parts, dest_ranks, comm, prev_reqs);
}

template <typename T>
void epa_mpi_recieve_merge(T& obj, std::vector<int>& src_ranks, MPI_Comm comm)
{
  for (auto rank : src_ranks)
  {
    T remote_obj;
    epa_mpi_recieve(remote_obj, rank, comm);
    merge(obj, remote_obj);
  }
}


template <typename T>
void epa_mpi_gather(T& obj, int dest_rank, std::vector<int>& src_ranks, int local_rank)
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
      epa_mpi_recieve(remote_obj, src_rank, MPI_COMM_WORLD);
      merge(obj, remote_obj);
    }
  }
}

template <typename T>
void epa_mpi_bcast(T& obj, int src_rank, std::vector<int>& dest_ranks, int local_rank)
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
    epa_mpi_recieve(obj, src_rank, MPI_COMM_WORLD);
  }
}

#endif

