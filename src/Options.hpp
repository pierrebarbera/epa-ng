#ifndef EPA_OPTIONS_H_
#define EPA_OPTIONS_H_

class Options {

public:
  Options ()
    : opt_insertion_branches(false), opt_model(false), opt_branches(false), support_threshold(0.01) {};

  Options (bool opt_insertion_branches, bool opt_model, bool opt_branches)
    : opt_insertion_branches(opt_insertion_branches), opt_model(opt_model), opt_branches(opt_branches)
    , support_threshold(0.01) {};

  ~Options () = default;

  bool opt_insertion_branches;
  bool opt_model;
  bool opt_branches;
  double support_threshold;

private:

};

#endif
