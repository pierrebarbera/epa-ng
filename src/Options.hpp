#pragma once

class Options {

public:
  Options ()
    : prescoring(false), opt_model(false), opt_branches(false), sliding_blo(true)
    , support_threshold(0.01), acc_threshold(false), prescoring_by_percentage(false)
    , prescoring_threshold(0.95), ranged(false), dump_binary_mode(false)
    , load_binary_mode(false), chunk_size(1000)
    { }

  ~Options () = default;

  bool prescoring;
  bool opt_model;
  bool opt_branches;
  bool sliding_blo;
  double support_threshold;
  bool acc_threshold;
  bool prescoring_by_percentage;
  double prescoring_threshold;
  bool ranged;
  bool dump_binary_mode;
  bool load_binary_mode;
  unsigned int chunk_size;
};
