#include "net/epa_mpi_util.hpp"

#include "net/mpihead.hpp"
#include "util/logging.hpp"

#include <limits>
#include <utility>
#include <cstddef>

std::pair<size_t, size_t> local_seq_package( const size_t num_seqs )
{
  int local_rank = 0;
  int num_ranks = 1;

  MPI_COMM_RANK(MPI_COMM_WORLD, &local_rank);
  MPI_COMM_SIZE(MPI_COMM_WORLD, &num_ranks);

  LOG_INFO << "Number of MPI ranks: " << num_ranks;

  size_t local_rank_seq_offset = 0;
  size_t part_size = num_seqs;

  if (num_ranks > 1) {
    // how many should each rank read?
    part_size = ceil(num_seqs / static_cast<double>(num_ranks));
    LOG_INFO << "Number of sequences per MPI rank: " << part_size;

    local_rank_seq_offset = part_size * local_rank;
  }
  return std::make_pair( local_rank_seq_offset, part_size );
}