#pragma once

#include "mpihead.hpp"
#include "streams.hpp"

#include <fstream>
#include <iostream>
#include <string>
#include <memory>
#include <streambuf>

class Log {
public:
  Log(const std::string& log_file, bool stdout=true) 
    : log_file_(new std::ofstream(log_file)),
      out_( new Teed_ostream(stdout ? std::cout.rdbuf() : nullptr, (*log_file_).rdbuf()) ),
      err_( new Prefixed_ostream(out_, "[ PEPA ERROR ] ") ),
      dbg_( new Prefixed_ostream(out_, "[ PEPA DEBUG ] ") )
    { };
  
  Log() : Log("/dev/null") { };

  Log(bool stdout) : Log("/dev/null", stdout) { };

  ~Log() 
  {
    if(log_file_) 
      log_file_->close();
  };

  void flush() {std::cout.flush();if(log_file_)(*log_file_).flush();};

  Log& operator=(Log&& other)
  {
    std::swap(log_file_, other.log_file_);
    std::swap(out_, other.out_);
    std::swap(err_, other.err_);
    std::swap(dbg_, other.dbg_);
    return *this;
  }

  template<typename T>
  Log& operator<<(T&& msg)
  {
    *out_ << msg;
    return *this;
  }

  Log& operator<<(std::ostream&(*f)(std::ostream&))
  {
    *out_ << f;
    return *this;
  }

  Prefixed_ostream& err() {return *err_;};
  Prefixed_ostream& dbg() {return *dbg_;};

private:
  std::unique_ptr<std::ofstream> log_file_ = nullptr;
  std::shared_ptr<Teed_ostream> out_;
  std::unique_ptr<Prefixed_ostream> err_;
  std::unique_ptr<Prefixed_ostream> dbg_;
};


extern Log lgr;
