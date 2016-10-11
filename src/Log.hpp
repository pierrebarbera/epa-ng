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
  Log(const std::string& log_file, bool stdout=true);
  Log();
  Log(bool stdout);
  ~Log();
  
  void flush();

  Log& operator=(Log&& other);
  
  template<typename T>
  Log& operator<<(T&& msg)
  {
    *out_ << msg;
    return *this;
  }

  Log& operator<<(std::ostream&(*f)(std::ostream&));
  
  Prefixed_ostream& err() {return *err_;}
  Prefixed_ostream& dbg() {return *dbg_;}

private:
  std::unique_ptr<std::ofstream> log_file_ = nullptr;
  std::shared_ptr<Teed_ostream> out_;
  std::unique_ptr<Prefixed_ostream> err_;
  std::unique_ptr<Prefixed_ostream> dbg_;
};


extern Log lgr;
