#pragma once

#ifdef __MPI
#include <mpi.h>
#define MPI_INIT(argc, argv) MPI_Init(argc, argv)
#define MPI_FINALIZE() MPI_Finalize()
#define MPI_COMM_RANK(comm, rank) MPI_Comm_rank(comm, rank)
#else
#define MPI_INIT(argc, argv)
#define MPI_FINALIZE()
#define MPI_COMM_RANK(comm, rank)
#endif
