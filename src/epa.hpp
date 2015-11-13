#include <string>
#include <vector>
#include <stdexcept>
#include <iostream>

#include "Tree.hpp"
#include "MSA.hpp"
#include "file_io.hpp"

void check_clvs_zero(pll_partition_t * p, pll_utree_t * t)
{
  // TODO garbage
  int tip_nodes_count = 8;
  int inner_nodes_count = 6;
  int i;
  pll_utree_t ** tipnodes = (pll_utree_t  **)calloc(tip_nodes_count,
                                                    sizeof(pll_utree_t *));
  pll_utree_query_tipnodes(t, tipnodes);

  pll_utree_t ** innernodes = (pll_utree_t  **)calloc(inner_nodes_count,
                                                    sizeof(pll_utree_t *));

  pll_utree_query_innernodes(t, innernodes);

  double logl = 0.0;
  for (i = 0; i <tip_nodes_count + inner_nodes_count; ++i)
  {
    if (i < tip_nodes_count)
    {
      t = tipnodes[i]->back;
    }
    else
    {
      t = innernodes[i - tip_nodes_count];
    }

    logl = pll_compute_edge_loglikelihood(p,
                                         t->clv_index,
                                         t->scaler_index,
                                         t->back->clv_index,
                                         t->back->scaler_index,
                                         t->pmatrix_index,
                                         0);
    std::cout << logl << std::endl;
  }
  // (void) t;
  // for (int i = 0; i < p->clv_buffers; ++i)
  // {
  //   for (int j = 0; j < p->sites; ++j) {
  //     if ((p->clv[i])[j] == 0.0)
  //       throw std::runtime_error{std::string("CLV found to be zero. Index: ") + std::to_string(i) +
  //                                             std::string(" ; ") + std::to_string(j)};
  //   }
  // }
};

static void epa(std::string& tree_file, std::string& reference_msa_file, std::string& query_msa_file,
                std::vector<double> base_frequencies, std::vector<double> substitution_rates)
{
	// sanitize input, detect file formats

	// call tree object, it builds itself from file
	auto model = Model(base_frequencies, substitution_rates);
	auto tree = Tree(tree_file, reference_msa_file, model);

  auto query_reads = build_MSA_from_file(query_msa_file);

  //tree.visit(&check_clvs_zero);

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
