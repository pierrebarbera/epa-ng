#include "epa.hpp"

#include <stdexcept>
#include <iostream>
#include <fstream>

#include "file_io.hpp"
#include "jplace_util.hpp"

using namespace std;

void epa(string& tree_file, string& reference_msa_file, string& query_msa_file,
                vector<double> base_frequencies, vector<double> substitution_rates,
                double alpha, string invocation)
{
	// sanitize input, detect file formats
  // TODO detect filetype

	// Build the reference tree
	auto model = Model(base_frequencies, substitution_rates, alpha);
  auto ref_msa = build_MSA_from_file(reference_msa_file);

  MSA query_msa;
  if (query_msa_file.size() != 0)
    query_msa = build_MSA_from_file(query_msa_file);

  auto tree = Tree(tree_file, ref_msa, model, query_msa);

  // place query sequences
  auto placements = tree.place();

  cout << placement_set_to_jplace_string(placements, invocation);

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
