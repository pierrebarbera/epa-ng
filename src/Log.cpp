#include "Log.hpp"

Log::Log(const std::string& log_file, bool stdout) 
  : log_file_(new std::ofstream(log_file)),
    out_( new Teed_ostream(stdout ? std::cout.rdbuf() : nullptr, (*log_file_).rdbuf()) ),
    err_( new Prefixed_ostream(*out_, "[ PEPA ERROR ] ") ),
    dbg_( new Prefixed_ostream(*out_, "[ PEPA DEBUG ] ") )
  { }

Log::Log() : Log("/dev/null") { }

Log::Log(bool stdout) : Log("/dev/null", stdout) { }

Log::~Log() 
{
  if(log_file_) 
    log_file_->close();
}

void Log::flush() {std::cout.flush();if(log_file_)(*log_file_).flush();};

Log& Log::operator=(Log&& other)
{
  std::swap(log_file_, other.log_file_);
  std::swap(out_, other.out_);
  std::swap(err_, other.err_);
  std::swap(dbg_, other.dbg_);
  return *this;
}

Log& Log::operator<<(std::ostream&(*f)(std::ostream&))
{
  *out_ << f;
  return *this;
}
