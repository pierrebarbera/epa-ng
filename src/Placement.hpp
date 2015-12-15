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
  inline double lwr() const {return lwr_;};
  inline double likelihood() const {return likelihood_;};
  inline double pendant_length() const {return pendant_length_;};
  inline double distal_length() const {return distal_length_;};
  inline unsigned int branch_id() const {return branch_id_;};

  // setters
  inline void lwr(double value) {lwr_ = value;};

private:
  unsigned int branch_id_;
  double likelihood_;
  double lwr_;
  double pendant_length_;
  double distal_length_;
};

#endif
