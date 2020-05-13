#include "sample/Placement.hpp"

Placement::Placement( Slim_Placement const& other )
    : branch_id_( other.branch_id() )
    , likelihood_( other.likelihood() )
    , lwr_( 0.0 )
    , pendant_length_( 0.0 )
    , distal_length_( 0.0 )
{
}

Slim_Placement::Slim_Placement( Placement const& other )
    : branch_id_( other.branch_id() )
    , likelihood_( other.likelihood() )
{
}
