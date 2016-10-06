#pragma once

#ifdef __MPI

#include "mpihead.hpp"
#include "Sample.hpp"

#include <sstream>
#include <memory>
#include <cereal/archives/binary.hpp>
#include <cereal/types/vector.hpp>
#include <cereal/types/string.hpp>

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
void epa_mpi_recieve(T& obj, int source_rank, MPI_Comm comm)
{
  // probe to find out the message size
  MPI_Status status;
  int size;
  err_check(MPI_Probe(source_rank, 0, comm, &status));
  MPI_Get_count(&status, MPI_CHAR, &size);

  // prepare buffer
  auto buffer = new char[size];

  //  get the actual payload
  err_check(MPI_Recv(buffer, size, MPI_CHAR, source_rank, 0, comm, &status));

  // deserialization
  std::stringstream ss;
  ss.write(buffer, size);
  cereal::BinaryInputArchive in_archive(ss);

  // build the object
  in_archive(obj);

  delete[] buffer;
}

template <typename T>
void epa_mpi_gather(T& obj, int dest_rank, std::unordered_map<int, int>& src_ranks, int local_rank)
{
  if (dest_rank != local_rank)
  {
    epa_mpi_send(obj, dest_rank, MPI_COMM_WORLD);
  }
  else
  {
    for (auto rank_pair : src_ranks)
    {
      int src_rank = rank_pair.second;
      
      if (local_rank == src_rank)
        continue;
      
      T remote_obj;
      epa_mpi_recieve(remote_obj, src_rank, MPI_COMM_WORLD);
      merge(obj, remote_obj);
    }
  }
}

template <typename T>
void epa_mpi_bcast(T& obj, int src_rank, std::unordered_map<int, int>& dest_ranks, int local_rank)
{
  if (src_rank == local_rank)
  {
    for (auto rank_pair : dest_ranks)
    {
      int dest_rank = rank_pair.second();
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

