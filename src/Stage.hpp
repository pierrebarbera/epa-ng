#pragma once

#include <functional>

#include "Token.hpp"

/**
 * Interface for Stages
 */
class Stage
{
// private:
//   Stage() = default;
//   virtual ~Stage() = default;

public:
  virtual void accept(Token&) = 0;
  virtual Token process(Token&) = 0;
  virtual void put(Token&) = 0;

  virtual bool exec() const final
  {
    return exec_;
  }

  virtual void exec(const bool b) final
  {
    exec_ = b;
  }

private:
  bool exec_ = true;
  
};

/**
 * Templated interface for Stages
 */
template <typename in_type, typename out_type>
class Typed_Stage : public Stage
{
public:
  using accept_func_type  = std::function<void(in_type)>;
  using process_func_type = std::function<out_type(in_type)>;
  using put_func_type     = std::function<void(out_type)>;

  Typed_Stage(accept_func_type   accept,
              process_func_type  process,
              put_func_type      put)
    : accept_(accept)
    , process_(process)
    , put_(put)
  { }

  Typed_Stage()		= delete;
  ~Typed_Stage() 	= default;

  void accept(Token& t) override
  {
    accept_(static_cast<in_type&>(t));
  }

  Token process(Token& token) override
  {
    return process_(static_cast<in_type&>(token));
  }

  void put(Token& t) override
  {
    put_(static_cast<out_type&>(t));
  }

  void set_accept_func(accept_func_type& f)
  {
    accept_(f);
  }

  void set_put_func(put_func_type& f)
  {
    put_(f);
  }
  

private:
  accept_func_type  accept_;
  process_func_type process_;
  put_func_type     put_;

};