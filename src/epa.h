#include "tree.h"

#include <string>

static void epa(std::string* tree_file, std::string* msa_file)
{
	(void) tree_file;
	(void) msa_file;
	// sanitize input, detect file formats
	
	// call tree object, it builds itself from file

	// for scalability... 
		// stream? then probably conjuntion between pll datatype and my wrapper needed
		// partition? further heuristic!

	// subtask handling!

	// for basic approach: call kernel directly, basic for loop, later expand to par solution
		// kernel does one thing: place one sequence on one branch and recalculate
		// do for every branch
			// do for every sequence
	// loop returns list of list of placements
	
	// filter: what fraction of placements to remember (absolute number, top x %, threshold)

	// pass placement lists to output writer, possibly also strategy based?
}
