#pragma once

#include <string>
#include <future>
#include <memory>
#include <sstream>

#include "sample/Sample.hpp"
#include "util/logging.hpp"
#include "io/jplace_util.hpp"

#ifdef __MPI
#include "net/epa_mpi_util.hpp"
#endif

class jplace_writer
{
public:
  jplace_writer()
  {
    init_mpi_();
  }
  jplace_writer(const std::string& out_dir,
                const std::string& file_name,
                const std::string& tree_string,
                const std::string& invocation_string)
  {
    init_mpi_();
    init_(out_dir, file_name, invocation_string);
    if (enclosed_) {
      *file_ << init_jplace_string(tree_string);
    }
  }

  ~jplace_writer()
  {
    // ensure last write/gather was completed
    wait();

    if (file_) {
      // finalize and close
      if (enclosed_) {
        *file_ << finalize_jplace_string(invocation_);
      }
      file_->close();
    }
  }

  jplace_writer& operator=( jplace_writer&& other )
  {
    this->invocation_ = std::move(other.invocation_);
    this->file_ = std::move(other.file_);
    other.file_ = nullptr;
    this->prev_gather_ = std::move(other.prev_gather_);
    return *this;
  }

  void write( Sample<>& chunk )
  {
    write_gather_(chunk);
  }

  void wait()
  {
    #ifdef __PREFETCH
    if (prev_gather_.valid()) {
      prev_gather_.get();
    }
    #endif
  }

protected:
  void write_gather_( Sample<>& chunk )
  {
    #ifdef __PREFETCH
      // ensure the last write has finished
      if (prev_gather_.valid()) {
        prev_gather_.get();
      }
      prev_gather_ = std::async(std::launch::async,
      [chunk = chunk, this]() mutable {
        this->gather_(chunk);
        this->write_(chunk);
      });
    #else
      gather_(chunk);
      write_(chunk);
    #endif
  }

  void write_(Sample<>& chunk)
  {
    if (file_) {
      if (first_){
        first_ = false;
      } else {
        *file_ << ",\n";
      }

      *file_ << sample_to_jplace_string(chunk);
    }
  }

  virtual void gather_( Sample<>& chunk )
  {
    #ifdef __MPI
    Timer<> dummy;
    epa_mpi_gather(chunk, 0, all_ranks_, local_rank_, dummy);
    #else
    (void) chunk;
    #endif //__MPI
  }

  void init_( const std::string& out_dir,
              const std::string& file_name,
              const std::string& invocation_string,
              const bool enclosed=true)
  {
    const auto file_path = out_dir + file_name;
    invocation_ = invocation_string;
    file_ = std::make_unique<std::fstream>();
    file_->open(file_path,
                std::fstream::in | std::fstream::out | std::fstream::trunc);

    if (not file_->is_open()) {
      throw std::runtime_error{file_name + ": could not open!"};
    }
    enclosed_ = enclosed;
  }

  void init_mpi_()
  {
    #ifdef __MPI // then have one outfile per rank
    int num_ranks = 0;
    MPI_COMM_RANK(MPI_COMM_WORLD, &local_rank_);
    MPI_COMM_SIZE(MPI_COMM_WORLD, &num_ranks);

    all_ranks_.resize(num_ranks);
    for (int i = 0; i < num_ranks; ++i) {
      all_ranks_[i] = i;
    }
    #endif
  }

  // void invocation(const std::string& invocation) {invocation_ = invocation;}
  // const std::string& invocation() const {return invocation_;}

protected:
  std::string invocation_;
  std::unique_ptr<std::fstream> file_ = nullptr;
  std::future<void> prev_gather_;
  bool first_ = true;
  bool enclosed_ = true;
  int local_rank_ = 0;
  std::vector<int> all_ranks_ = {0};
};

/**
 * special jplace writer that writes part files per rank, then does a parallel read/write to merge them in the end
 */
class localized_jplace_writer : public jplace_writer
{
public:
  localized_jplace_writer() = default;
  localized_jplace_writer(const std::string& out_dir,
                          const std::string& file_name,
                          const std::string& tree_string,
                          const std::string& invocation_string)
  {
    auto local_file = file_name;
    #ifdef __MPI // then have one outfile per rank
    local_file = std::to_string(local_rank_) + "." + local_file;
    #endif

    init_(out_dir, local_file, invocation_string, false);

    final_file_ = out_dir + file_name;
    tree_string_ = tree_string;
  }


  ~localized_jplace_writer()
  {
    // ensure everything is dandy, then merge the files
    wait();

    #ifdef __MPI
    const bool master = (local_rank_ == 0);
    // size of this part
    size_t num_bytes = file_->tellp();

    // measure the size of the leading and trailing string
    std::string leading;
    std::string trailing;
    if (master) {
      leading   = init_jplace_string(tree_string_);
      trailing  = finalize_jplace_string(invocation_);
      // rank 0 writes the initial string into its fraction
      num_bytes += leading.size();
    }

    // broadcast all block sizes, such that everyone can calculate their displacement
    std::vector<size_t> block_sizes( all_ranks_.size() );
    MPI_Allgather(&num_bytes, 1, MPI_SIZE_T, block_sizes.data(), 1, MPI_SIZE_T, MPI_COMM_WORLD);
    size_t displacement = 0;
    for (int i = 0; i < local_rank_; ++i) {
      displacement += block_sizes[i];
    }

    // for (auto& )
    LOG_DBG << "Displacement: " << displacement;

    // create parallel outfile, with num_bytes reserved for this rank 
    MPI_File fh;
    MPI_Status status;
    MPI_File_open(MPI_COMM_WORLD, final_file_.c_str(), MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL,&fh);
    // MPI_File_set_view(fh, displacement, MPI_CHAR, filetype, datarep, MPI_INFO_NULL);
    
    // get the contents of the file back
    std::stringstream buffer;
    buffer << leading;
    file_->seekg(0); // back to the start
    buffer << file_->rdbuf();

    MPI_File_write_at_all(fh, displacement, buffer.str().c_str(), buffer.str().size(),
                              MPI_CHAR, &status);

    if (master) {
      // goto end of file
      size_t eof = std::accumulate(std::begin(block_sizes), std::end(block_sizes), 0);
      MPI_File_seek(fh, eof, MPI_SEEK_SET);
      MPI_File_write(fh, trailing.c_str(), trailing.size(),
                   MPI_CHAR, &status);
    }

    MPI_File_close(&fh);
    #endif
  }

private:
  // no gather, only write
  void gather_( Sample<>& chunk ) override
  {
    (void) chunk;
  }

  std::string final_file_;
  std::string tree_string_;

};
 