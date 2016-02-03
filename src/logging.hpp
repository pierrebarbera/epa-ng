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

  template<typename T>
  Log& operator<<(T&& msg)
  {
    std::cout << msg;
    if(log_file_)
      (*log_file_) << msg;
    return *this;
  }

  Log& operator<<(std::ostream&(*f)(std::ostream&))
  {
    std::cout << f;
    if(log_file_)
      (*log_file_) << f;
    return *this;
  }

private:
  std::unique_ptr<std::ofstream> log_file_ = nullptr;
};

extern Log lgr;
