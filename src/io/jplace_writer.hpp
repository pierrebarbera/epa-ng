#pragma once

#include <string>
#include <future>

#include "sample/Sample.hpp"
#include "util/logging.hpp"

#ifdef __MPI
#include "net/epa_mpi_util.hpp"
#endif

class Jplace_writer
{
public:
  Jplace_writer() = default;
  Jplace_writer(const std::string& file_name,
                const std::string& tree_string,
                const std::string& invocation_string)
    : invocation_(invocation_string)
  {
    file_.open(file_name);
    file_ << init_jplace_string(tree_string);
  }

  ~Jplace_writer()
  {
    if (file_) {
      // ensure last write/gather was completed
      #ifdef __PREFETCH
        if (prev_gather_.valid()) {
          prev_gather_.get();
        }
      #endif

      // finalize and close
      file_ << finalize_jplace_string(invocation_);
      file_.close();
    }
  }

  Jplace_writer& operator=(Jplace_writer&& other)
  {
    this->invocation_ = std::move(other.invocation_);
    this->file_ = std::move(other.file_);
    this->prev_gather_ = std::move(other.prev_gather_);
    return *this;
  }

  void gather_write(Sample<>& chunk,
                    const std::vector<int>& all_ranks={0},
                    const int local_rank=0)
  {
    if (file_) {
      #ifdef __PREFETCH
        // ensure the last write has finished
        if (prev_gather_.valid()) {
          prev_gather_.get();
        }
        prev_gather_ = std::async(std::launch::async,
        [chunk, all_ranks, local_rank, this]() mutable -> void{
          gather_(chunk, all_ranks, local_rank);
          write_(chunk);
        });
      #else
        gather_(chunk, all_ranks, local_rank);
        write_(chunk);
      #endif
    }
  }

private:
  void write_( Sample<>& chunk)
  {
    if (file_) {
      file_ << sample_to_jplace_string(chunk) << ",\n";
    }
  }

  void gather_( Sample<>& chunk,
                    const std::vector<int>& all_ranks,
                    const int local_rank)
  {
    #ifdef __MPI
    Timer<> dummy;
    epa_mpi_gather(chunk, 0, all_ranks, local_rank, dummy);
    #endif //__MPI
  }

  std::string invocation_;
  std::ofstream file_;
  std::future<void> prev_gather_;
};
 