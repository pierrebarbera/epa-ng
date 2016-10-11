#pragma once

#include <iostream>
#include <string>
#include <streambuf>
#include <memory>

class Teed_ostream 
{
public:
  Teed_ostream(std::streambuf * sb1, std::streambuf * sb2)
    : stream_a(sb1), stream_b(sb2)
    { };
  
  Teed_ostream() 
    : stream_a(nullptr), stream_b(nullptr)
    { };

  ~Teed_ostream() = default;

  template<typename T>
  Teed_ostream& operator<<(T&& msg)
  {
    stream_a << msg;
    stream_b << msg;
    return *this;
  }

  Teed_ostream& operator<<(std::ostream&(*f)(std::ostream&))
  {
    stream_a << f;
    stream_b << f;
    return *this;
  }
private:
  std::ostream stream_a;
  std::ostream stream_b;
};

class Prefixed_ostream
{
public:
  Prefixed_ostream(std::shared_ptr<Teed_ostream> sb, std::string prefix) 
    : stream_(sb), prefix_(prefix)
    { };

  Prefixed_ostream(Prefixed_ostream&& other)
  { 
    std::swap(stream_, other.stream_);
    std::swap(prefix_, other.prefix_);
  }

  Prefixed_ostream() : stream_(nullptr) { };

  ~Prefixed_ostream() = default;

  Prefixed_ostream& operator=(Prefixed_ostream&& other)
  {
    std::swap(stream_, other.stream_);
    std::swap(prefix_, other.prefix_);
    return *this;
  }

  template<typename T>
  Prefixed_ostream& operator<<(T&& msg)
  {
    *stream_ << prefix_ << msg << std::endl;
    return *this;
  }

  Prefixed_ostream& operator<<(std::ostream&(*f)(std::ostream&))
  {
    *stream_ << prefix_ << f << std::endl;
    return *this;
  }

private:
  std::shared_ptr<Teed_ostream> stream_;
  std::string prefix_;
};
