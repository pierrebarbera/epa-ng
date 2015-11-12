#include <string>
#include <vector>

#include "Tree.h"
#include "MSA.h"
#include "file_io.h"

static void epa(std::string& tree_file, std::string& reference_msa_file, std::string& query_msa_file,
                std::vector<double> base_frequencies, std::vector<double> substitution_rates)
{
	// sanitize input, detect file formats

	// call tree object, it builds itself from file
	auto model = Model(base_frequencies, substitution_rates);
	auto tree = Tree(tree_file, reference_msa_file, model);

  auto query_reads = build_MSA_from_file(query_msa_file);


  tree.place(query_reads.get(0));

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
