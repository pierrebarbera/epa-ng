#ifndef EPA_PLACEMENT_H_
#define EPA_PLACEMENT_H_

class Placement {
public:
  Placement() = delete;

  Placement(unsigned int branch_id, double likelihood, double pendant_length,
    double distal_length)
    : branch_id_(branch_id), likelihood_(likelihood), lwr_(0.0), pendant_length_(pendant_length),
      distal_length_(distal_length)
  {};

  Placement(Placement const& other) = default;
  Placement(Placement&& other) = default;

  Placement& operator = (Placement const& other) = default;
  Placement& operator = (Placement && other) = default;

  ~Placement() = default;

  // getters
  double lwr() const {return lwr_;};
  double likelihood() const {return likelihood_;};
  double pendant_length() const {return pendant_length_;};
  double distal_length() const {return distal_length_;};
  unsigned int branch_id() const {return branch_id_;};

  // setters
  void lwr(double value) {lwr_ = value;};
  void likelihood(double value) {likelihood_ = value;};
  void pendant_length(double value) {pendant_length_ = value;};
  void distal_length(double value) {distal_length_ = value;};

private:
  unsigned int branch_id_;
  double likelihood_;
  double lwr_;
  double pendant_length_;
  double distal_length_;
};

#endif
