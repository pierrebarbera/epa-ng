#pragma once

#include <fstream>
#include <iostream>
#include <string>
#include <memory>

class Log {
public:
  Log(const std::string& log_file) : log_file_(new std::ofstream(log_file)) {};
  Log() = default;
  ~Log() {if(log_file_) log_file_->close();};

  Log& operator=(Log&& other)
  {
    std::swap(log_file_, other.log_file_);
    return *this;
  }

  // TODO templetize?
  std::ostream& operator<<(std::string&& msg)
  {
    std::cout << msg << std::endl;
    if(log_file_)
      (*log_file_) << msg << std::endl;
    return std::cout;
  }

private:
  std::unique_ptr<std::ofstream> log_file_ = nullptr;
};

Log lgr;
