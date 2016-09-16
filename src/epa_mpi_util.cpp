#include "epa_mpi_util.hpp"

#include <sstream>
#include <memory>
#include <cereal/archives/binary.hpp>

#ifdef __MPI
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

void epa_mpi_send(Sample& sample, int dest_rank, MPI_Comm comm)
{
  // serialize the Sample
  std::stringstream ss;
  cereal::BinaryOutputArchive out_archive(ss);
  out_archive(sample);

  // send sample to specified node
  std::string data = ss.str();
  std::unique_ptr<char> buffer(new char[data.size()]);
  memcpy(buffer.get(), data.c_str(), data.size() * sizeof(char));
  err_check(MPI_Send(buffer.get(), data.size(), MPI_CHAR, dest_rank, 0, comm));

}

void epa_mpi_recieve(Sample& sample, int source_rank, MPI_Comm comm)
{
  // probe to find out the message size
  MPI_Status status;
  int size;
  err_check(MPI_Probe(source_rank, 0, comm, &status));
  MPI_Get_count(&status, MPI_CHAR, &size);

  // prepare buffer
  std::unique_ptr<char> buffer(new char[size]);

  //  get the actual payload
  err_check(MPI_Recv(buffer.get(), size, MPI_CHAR, source_rank, 0, comm, &status));

  // deserialization
  std::stringstream ss;
  ss.write(buffer.get(), size);
  cereal::BinaryInputArchive in_archive(ss);

  // build the sample object
  in_archive(sample);
}
#endif
