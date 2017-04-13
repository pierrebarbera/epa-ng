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
#include "function_traits.hpp"
#include "template_magic.hpp"

/**
 * Building a Stage Tuple out of a bunch of lambda functions/functors
 */
template < class I, class... lambdas>
struct stage_types_base;

template < std::size_t... I, class... lambdas >
struct stage_types_base<std::index_sequence<I...>, lambdas...>
{
  using types = typename std::tuple< Typed_Stage<I, lambdas>... >;
};

template < class... lambdas >
struct stage_types 
  : stage_types_base<std::make_index_sequence<sizeof...(lambdas) >, lambdas...>
{
};

/**
 * Basic Pipeline Class. Runs all stages in serial.
 */
template <class... lambdas>
class Pipeline
{
  using stack_type      = typename stage_types< lambdas... >::types;
  using token_set_type  = typename token_types< stack_type >::types;

public:
  using hook_type       = std::function<void()>;

  Pipeline(const stack_type & stages, const hook_type & per_loop_hook)
    : stages_(stages)
    , per_loop_hook_(per_loop_hook)
  { }

  ~Pipeline() = default;

  template <class Function>
  auto push(const Function& f) const
  {
    constexpr size_t num_stages = sizeof...(lambdas) + 1u;
    constexpr auto new_stage_id = num_stages - 1u;

    using stage_type = Typed_Stage<new_stage_id, Function>;
    using new_stack_type = typename stage_types<lambdas..., Function>::types;

    new_stack_type stage_tuple
      = std::tuple_cat(stages_, std::make_tuple(stage_type(f)));
    
    return Pipeline<lambdas..., Function>(stage_tuple, per_loop_hook_);
  }

  void process()
  {
    token_set_type tokens;

    // "last" token that is still used on the particular MPI-Rank (or thread or...)
    Token const * last_token = nullptr;

    do { 

      // per-loop pre-hook
      per_loop_hook_();

      for_each(stages_, [&](auto& s) {

        if (s.exec()) {

          constexpr auto stage_id = s.id();

          auto& in_token = std::get<stage_id>(tokens);
          auto& out_token = std::get<stage_id+1u>(tokens);
          
          s.accept(in_token); // noop if shared mem, mpi_merge_receive if mpi

          out_token = s.process(in_token); // casts asbtract token to its input type and runs the user code

          if (stage_id != 0) {
            // carry over the token status
            out_token.status(in_token.status());
          }

          s.put(out_token); // noop if shared mem, mpi_split_send if mpi

          last_token = &out_token;

        }
      });

      if(last_token->rebalance()) {
        // do the kansas city shuffle...
      }

    } while (last_token->valid()); //returns valid if data token or default initialized
  }

private:
  stack_type stages_;
  hook_type per_loop_hook_;

};

template <class stage_f>
auto make_pipeline( const stage_f& first_stage, 
                    const typename Pipeline<stage_f>::hook_type& per_loop_hook
                      = [](){}) 
{
  return Pipeline<stage_f>(std::make_tuple(Typed_Stage<0u, stage_f>(first_stage)), per_loop_hook);
}

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
