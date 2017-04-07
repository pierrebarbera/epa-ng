#pragma once

#include <functional>
#include <type_traits>

#include "Token.hpp"
#include "function_traits.hpp"


/**
 * Interface for Stages
 */
class Stage
{
// private:
//   Stage() = default;
//   virtual ~Stage() = default;

public:
  // virtual void accept(Token&) = 0;
  // virtual Token process(Token&) = 0;
  // virtual void put(Token&) = 0;

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
template <typename lambda>
class Typed_Stage : public Stage
{
public:
  using traits = thrill::common::FunctionTraits<lambda>;
  using base_in_type = typename traits::template arg<0>;
  using base_out_type = typename traits::result_type;
  using in_type = typename std::remove_reference<base_in_type>::type;
  using out_type = typename std::remove_reference<base_out_type>::type;
  using accept_func_type  = std::function<void(in_type&)>;
  using process_func_type = std::function<out_type(in_type&)>;
  using put_func_type     = std::function<void(out_type&)>;

  Typed_Stage(accept_func_type   accept,
              process_func_type  process,
              put_func_type      put)
    : accept_(accept)
    , process_(process)
    , put_(put)
  { }

  Typed_Stage(process_func_type process)
    : accept_([](in_type&) { })
    , process_(process)
    , put_([](out_type&) { })
  { }

  Typed_Stage()		= delete;
  ~Typed_Stage() 	= default;

  void accept(Token& arg)
  {
    accept_(static_cast<in_type&>(arg));
  }

  Token process(Token& arg)
  {
    return static_cast<Token>(process_(static_cast<in_type&>(arg)));
  }

  void put(Token& arg) 
  {
    put_(static_cast<out_type&>(arg));
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
