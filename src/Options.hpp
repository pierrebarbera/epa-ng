#pragma once

class Options {

public:
  Options ()
    : prescoring(false), opt_model(false), opt_branches(false), support_threshold(0.01)
    , acc_threshold(false), prescoring_by_percentage(false), prescoring_threshold(0.95) {};

  Options (bool prescoring, bool opt_model, bool opt_branches)
    : prescoring(prescoring), opt_model(opt_model), opt_branches(opt_branches)
    , support_threshold(0.01), acc_threshold(false), prescoring_by_percentage(false), prescoring_threshold(0.95) {};

  ~Options () = default;

  bool prescoring;
  bool opt_model;
  bool opt_branches;
  double support_threshold;
  bool acc_threshold;
  bool prescoring_by_percentage;
  double prescoring_threshold;

private:

};
