#pragma once

#include "epa_mpi_util.hpp"

template <class T>
class Communicator {
public:
  Communicator () = default;
  ~Communicator() = default;
public:
  virtual T get() = 0;
  virtual void put(T data) = 0;
};

template<class data_type>
class MPI_Communicator : public Communicator<data_type>  {
private:
  int source_rank;
  int dest_rank;
public:
  MPI_Communicator(int source_rank, int dest_rank) : source_rank(source_rank), dest_rank(dest_rank) {};
  ~MPI_Communicator() = default;

  data_type get() override
  {
    data_type remote_data;
    epa_mpi_recieve(remote_data, source_rank, MPI_COMM_WORLD);
    return remote_data;
  }

  void put(data_type data) override
  {
    epa_mpi_send(data, dest_rank, MPI_COMM_WORLD);
  }
};

template<class data_type>
class Memory_Communicator : public Communicator<data_type> {
private:
  data_type location;
public:
  Memory_Communicator(data_type location) : location(location) {};
  ~Memory_Communicator() = default;

  // for threaded version: mutexes?

  data_type get() override
  {
    return location;
  }

  void put(data_type data) override
  {
    location = data;
  }
};
