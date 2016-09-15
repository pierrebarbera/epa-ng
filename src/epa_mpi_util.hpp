#pragma once

#include "mpihead.hpp"
#include "Sample.hpp"

#ifdef __MPI
void epa_mpi_send(Sample& sample, int dest_rank, MPI_Comm comm);
void epa_mpi_recieve(Sample& sample, int source_rank, MPI_Comm comm);
#endif
