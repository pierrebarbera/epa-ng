#include "PQuery_Set.hpp"

using namespace std;

PQuery_Set::PQuery_Set()
{ }

PQuery_Set::PQuery_Set(const unsigned int size) : pquerys_(size)
{ }

PQuery_Set::PQuery_Set(const string newick) : newick_(newick)
{ }

PQuery_Set::~PQuery_Set ()
{ }
