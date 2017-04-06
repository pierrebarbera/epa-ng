#pragma once

#include <functional>
#include <vector>
#include <stdexcept>
#include <type_traits>
#include <tuple>
#include <memory>

#include "Stage.hpp"
#include "Token.hpp"
#include "schedule.hpp"
#include "function_stack.hpp"
#include "function_traits.hpp"
#include "template_magic.hpp"

/**
 * Basic Pipeline Class. Runs all stages in serial.
 */
// template <typename input, typename output>
class Pipeline
{
  using global_in  = std::function<void(Token&)>;
  using global_out = std::function<void(Token&)>;
  // static constexpr auto size_ = sizeof...(lambdas);

  static void noop(Token&) { }

public:
  template <class... lambdas>
  Pipeline( lambdas... fn )
    // : stages_(source_function, fn...)
  {
    constexpr size_t num_tasks = sizeof...(lambdas);
    static_assert(num_tasks >= 3, "Pipeline must be assigned at least 3 function types.");

    std::tuple<lambdas...> f_tuple{fn...};
    // in_func_ = head(f_tuple);
    // // just the innermost functions: skip the first one, then of that set skip the last one
    // auto process_funcs  = rskip(skip(f_tuple));
    // out_func_ = tail(f_tuple); //std::get<sizeof...(lambdas) - 1u>(f_tuple);

    size_t task_num = 0;
    for_each(f_tuple, [&](auto f) {
      using traits = thrill::common::FunctionTraits<decltype(f)>;
      using in_type = typename traits::template arg<0>;
      using out_type = typename traits::result_type;

      std::function<void(in_type)> in_func = [](in_type){};
      std::function<void(out_type)> out_func = [](out_type){};
      // global_out out_func = Pipeline::noop;
      // if ( task_num == 0 ) {
      //   in_func = in_func_;
      // } 
      // if ( task_num == num_tasks - 1u) {
      //   out_func = out_func_;
      // }
      stages_.emplace_back(
        new Typed_Stage<in_type, out_type>(
          in_func,
          f,
          out_func
        )
      );
      ++task_num;
    });

    for (auto& s : stages_) {
      s->exec(true);
    }
  }

  ~Pipeline() = default;

  void process()
  {
    Token token; 
    while (token) { //returns valid if data token or default initialized

      if(token.rebalance()) {
        // do the kansas city shuffle...
      }

      for (auto& s : stages_) {

        if (s->exec()) {
          
          s->accept(token); // noop if shared mem, mpi_merge_receive if mpi
          
          token = s->process(token); // casts asbtract token to its input type and runs the user code
          
          s->put(token); // noop if shared mem, mpi_split_send if mpi

        }
      }
    }
  }

private:
  std::vector<std::unique_ptr<Stage>> stages_;

};

#ifdef __MPI

#include "epa_mpi_util.hpp"

/**
 * MPI Pipeline Class. On each rank, runs the parts of the pipeline assigned to that rank
 */
template <typename input, typename output, typename ... lambdas>
class MPI_Pipeline : public Pipeline<input, output>
{
public:
  using super_type  = Pipeline<input, output>; 
  using source_type = super_type::source_type;
  using sink_type   = super_type::sink_type;

  MPI_Pipeline( std::vector<double> initial_difficulty,
                lambdas... fn)
    : super_type(fn...)
  {
    if (stages.size() != initial_difficulty.size()) {
      throw std::runtime_error{"Initial difficulty must match number of stages!"};
    }

    int local_rank = 0;
    int local_stage;
    int world_size;
    MPI_Comm_rank(MPI_COMM_WORLD, &local_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    // get initial schedule
    auto init_nps = solve(num_stages, world_size, initial_difficulty);
    assign(local_rank, init_nps, schedule_, &local_stage);

    lgr.dbg() << "Schedule: ";
    for (size_t i = 0; i < schedule.size(); ++i) {
      lgr.dbg() << schedule[i].size() << " ";
    }
    lgr.dbg() << std::endl;

    for (size_t stage = 0; stage < stages_.size(); ++stage) {
      stages_[stage].exec = (local_stage == stage);
    }


  }

  MPI_Pipeline() = delete;

  ~MPI_Pipeline() = default;

private:
  schedule_type schedule_;
};

#endif
