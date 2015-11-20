#include "epa.hpp"

#include <stdexcept>

#include "file_io.hpp"

using namespace std;

void epa(string& tree_file, string& reference_msa_file, string& query_msa_file,
                vector<double> base_frequencies, vector<double> substitution_rates,
                double alpha)
{
	// sanitize input, detect file formats

	// call tree object, it builds itself from file
	auto model = Model(base_frequencies, substitution_rates, alpha);
	auto tree = Tree(tree_file, reference_msa_file, model);

  auto query_reads = build_MSA_from_file(query_msa_file);

  tree.place(query_reads);

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
