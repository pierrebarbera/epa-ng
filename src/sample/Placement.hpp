#pragma once

#include <cstddef>

class Preplacement;

class Placement {
  public:
  using branch_id_t = unsigned int;

  Placement() = default;

  Placement( branch_id_t const branch_id,
             double const likelihood,
             double const pendant_length,
             double const distal_length )
      : branch_id_( branch_id )
      , likelihood_( likelihood )
      , lwr_( 0.0 )
      , pendant_length_( pendant_length )
      , distal_length_( distal_length )
  {
  }

  Placement( Preplacement const& other );

  Placement( Placement const& other ) = default;
  Placement( Placement&& other )      = default;

  Placement& operator=( Placement const& other ) = default;
  Placement& operator=( Placement&& other ) = default;

  ~Placement() = default;

  // getters
  double lwr() const { return lwr_; }
  double likelihood() const { return likelihood_; }
  double pendant_length() const { return pendant_length_; }
  double distal_length() const { return distal_length_; }
  branch_id_t branch_id() const { return branch_id_; }

  // setters
  void lwr( double value ) { lwr_ = value; }
  void likelihood( double value ) { likelihood_ = value; }
  void pendant_length( double value ) { pendant_length_ = value; }
  void distal_length( double value ) { distal_length_ = value; }

  // serialization
  template< class Archive > void serialize( Archive& ar )
  {
    ar( branch_id_, likelihood_, lwr_, pendant_length_, distal_length_ );
  }

  private:
  branch_id_t branch_id_;
  double likelihood_;
  double lwr_;
  double pendant_length_;
  double distal_length_;
};

class Preplacement {
  public:
  using branch_id_t = unsigned int;
  Preplacement()    = default;

  Preplacement( branch_id_t const branch_id, double const likelihood )
      : branch_id_( branch_id )
      , likelihood_( likelihood )
      , lwr_( 0.0 )
  {
  }
  Preplacement( Placement const& other );

  Preplacement( Preplacement const& other ) = default;
  Preplacement( Preplacement && other ) = default;

  ~Preplacement() = default;

  Preplacement& operator=( Preplacement const& other ) = default;
  Preplacement& operator=( Preplacement&& other ) = default;

  // getters
  branch_id_t branch_id() const { return branch_id_; }
  double likelihood() const { return likelihood_; }
  double lwr() const { return lwr_; }

  // setters
  void likelihood( double value ) { likelihood_ = value; }
  void lwr( double value ) { lwr_ = value; }

  // serialization
  template< class Archive > void serialize( Archive& ar )
  {
    ar( branch_id_, likelihood_ );
  }

  private:
  branch_id_t branch_id_;
  double likelihood_;
  double lwr_;
};
