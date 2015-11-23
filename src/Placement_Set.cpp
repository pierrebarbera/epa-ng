#include "Placement_Set.hpp"

using namespace std;

Placement_Set::Placement_Set()
{ }

Placement_Set::Placement_Set(const unsigned int size) : placements_(size)
{ }

Placement_Set::Placement_Set(const string newick) : newick_(newick)
{ }

Placement_Set::~Placement_Set ()
{ }
