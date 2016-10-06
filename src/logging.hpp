#pragma once

#include "mpihead.hpp"

#include <fstream>
#include <iostream>
#include <string>
#include <memory>

class Log {
public:
  Log(const std::string& log_file) : log_file_(new std::ofstream(log_file)) {};
  Log() = default;
  ~Log() {if(log_file_) log_file_->close();};

  void to_cout(bool b) {to_cout_ = b;};
  void flush() {std::cout.flush();(*log_file_).flush();};

  Log& operator=(Log&& other)
  {
    std::swap(log_file_, other.log_file_);
    return *this;
  }

  template<typename T>
  Log& operator<<(T&& msg)
  {
    #ifdef __MPI
    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank==0)
    #endif
    {
    if (to_cout_)
      std::cout << msg;
    }
    if(log_file_)
      (*log_file_) << msg;
    return *this;
  }

  Log& operator<<(std::ostream&(*f)(std::ostream&))
  {
    #ifdef __MPI
    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank==0)
    #endif
    {
    if (to_cout_)
      std::cout << f;
    }
    if(log_file_)
      (*log_file_) << f;
    return *this;
  }

private:
  std::unique_ptr<std::ofstream> log_file_ = nullptr;
  bool to_cout_ = true;
};

extern Log lgr;
