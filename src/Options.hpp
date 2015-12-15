#ifndef EPA_OPTIONS_H_
#define EPA_OPTIONS_H_

class Options {
public:
  Options () : opt_insertion_branches(false), opt_model(false), opt_branches(false) {};
  Options (bool opt_insertion_branches, bool opt_model, bool opt_branches)
  : opt_insertion_branches(opt_insertion_branches), opt_model(opt_model), opt_branches(opt_branches) {};
  ~Options () = default;

  bool opt_insertion_branches;
  bool opt_model;
  bool opt_branches;

private:

};

#endif
