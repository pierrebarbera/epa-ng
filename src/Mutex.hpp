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

class Scoped_Mutex
{
public:
  explicit Scoped_Mutex(Mutex& m) 
    : mutex_(m)
  { 
    mutex_.lock(); 
  }

  ~Scoped_Mutex() 
  { 
    mutex_.unlock(); 
  }

private:
 Mutex& mutex_;
 bool locked_;

 void operator=(const Scoped_Mutex&);
 Scoped_Mutex(const Scoped_Mutex&);
};

class Mutex_List
{
public:
  Mutex_List(const size_t num_locks) : locks_(num_locks) { }
  Mutex_List() = default;
  ~Mutex_List() = default;

  Mutex& operator[] (const size_t index) { return locks_[index]; }
private:
  std::vector<Mutex> locks_;
  
};
