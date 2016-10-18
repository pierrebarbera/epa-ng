#pragma once

#ifdef __OMP
# include <omp.h>
class Mutex
{
public:
 Mutex() { omp_init_lock(&lock_); }
 ~Mutex() { omp_destroy_lock(&lock_); }
 void lock() { omp_set_lock(&lock_); }
 void unlock() { omp_unset_lock(&lock_); }
 
 Mutex(const Mutex& ) { omp_init_lock(&lock_); }
 Mutex& operator= (const Mutex& ) { return *this; }
private:
 omp_lock_t lock_;
};
#else

class Mutex
{
public:
 void lock() {}
 void unlock() {}
};
#endif

class Named_Lock
{
public:
  Named_Lock(const size_t num_locks) : locks_(num_locks) { }
  Named_Lock() = default;
  ~Named_Lock() = default;

  Mutex& operator[] (const size_t index) { return locks_[index]; }
private:
  std::vector<Mutex> locks_;
  
};
