#pragma once

#include <memory>
#include <vector>

#ifdef __OMP
#include <omp.h>

class Mutex
{
public:
 Mutex() { omp_init_lock(&lock_); }
 ~Mutex() { omp_destroy_lock(&lock_); }
 
 Mutex(const Mutex& ) { omp_init_lock(&lock_); }
 Mutex& operator= (const Mutex& ) { return *this; }

 void lock() { omp_set_lock(&lock_); }
 void unlock() { omp_unset_lock(&lock_); }
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
  explicit Scoped_Mutex(std::shared_ptr<Mutex>& m) : mutex_(m) { mutex_->lock(); }
  ~Scoped_Mutex() { mutex_->unlock(); }
  
  Scoped_Mutex(const Scoped_Mutex&) = delete;
  void operator=(const Scoped_Mutex&) = delete;
private:
	std::shared_ptr<Mutex> mutex_;
};

class Mutex_List
{
public:
  Mutex_List(const size_t num_locks)
  {
  	for (size_t i = 0; i < num_locks; ++i)
  	{
  		locks_.emplace_back(new Mutex);
  	}
  }

  Mutex_List() = default;
  ~Mutex_List() = default;

  std::shared_ptr<Mutex>& operator[] (const size_t index) { return locks_[index]; }
private:
  std::vector<std::shared_ptr<Mutex>> locks_;
};
