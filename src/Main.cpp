#include <iostream>
#include <string>
#include "epa.h"

#ifndef __PLL__
#define __PLL__
#include "pll.h"
#endif

using namespace std;

void plltest();

int main(int argc, char** argv)
{
	plltest();

	(void) argc;
	(void) argv;
	// parse inputs

	// epa call: object generation (static functions in epa.h?), includes sanitation
	auto tree_file = new string("abc");
	auto msa_file = new string("def");
	epa(tree_file, msa_file);
	// status callback?
	printf("gday\n");
 return 0;
}

void plltest()
{
	pll_partition_t * tree;

	tree = pll_partition_create(
		4, //num tips
		3, //num extra clv buf = inner nodes
		4, //num states = 4 for DNA
		6, //sequence length
		1, //num different subst. models
		3, //num prob matrices (2xleaves - 2), bzw num unique branch lengths?
		4, //num rate categories
		3, //num scale buffers = inner nodes
		PLL_ATTRIB_ARCH_SSE);

	/* initialize an array of frequencies */
 double frequencies[4] = { 0.25, 0.25, 0.25, 0.25 };
 /* set frequencies */
 pll_set_frequencies(tree, 0, frequencies);

	/* substitution rates for the GTR model */
 double subst_params[6] = {1,1,1,1,1,1}; //6 because (states*(states-1))/2
 /* set substitution parameters */
 pll_set_subst_params(tree, 0, subst_params);

 double rate_cats[4] = {	0.13695378267140107,  
                         0.47675185617665189,  
                         0.99999999997958422,  
                         2.38629436117236260};
 /* set rate categories */
 pll_set_category_rates(tree, rate_cats);

 /* setting the tip "states", CLVs, but specified as sequences with appropriate conversion function map_nt */
 pll_set_tip_states(tree, 0, pll_map_nt, "CACACA");
 pll_set_tip_states(tree, 1, pll_map_nt, "CACACD");
 pll_set_tip_states(tree, 2, pll_map_nt, "AGGACA");
 pll_set_tip_states(tree, 3, pll_map_nt, "CGTAGT");

 /* set unique branch lengths and their future indices f prob. matrices */
 double branch_lengths[3] = { 0.4, 0.5, 0.6 };
 int matrix_indices[3] = { 0, 1, 2};

 pll_update_prob_matrices(tree, 0, matrix_indices, branch_lengths, 3);

 /* operations structure, used to tell how to traverse the tree (inner node config)*/
 pll_operation_t * ops;

 ops = (pll_operation_t *)malloc(3 * sizeof(pll_operation_t));

 ops[0].parent_clv_index    = 4; //new index
 ops[0].parent_scaler_index = 0; //using which scaler
 ops[0].child1_clv_index    = 0; //child indexes
 ops[0].child2_clv_index    = 1;
 ops[0].child1_matrix_index = 1; //branch length prom matrix index (see branch_lengths)
 ops[0].child2_matrix_index = 1;
 ops[0].child1_scaler_index = PLL_SCALE_BUFFER_NONE;
 ops[0].child2_scaler_index = PLL_SCALE_BUFFER_NONE;

 ops[1].parent_clv_index    = 5;
 ops[1].parent_scaler_index = 1;
 ops[1].child1_clv_index    = 2;
 ops[1].child2_clv_index    = 3;
 ops[1].child1_matrix_index = 2;
 ops[1].child2_matrix_index = 2;
 ops[1].child1_scaler_index = PLL_SCALE_BUFFER_NONE;
 ops[1].child2_scaler_index = PLL_SCALE_BUFFER_NONE;

	ops[2].parent_clv_index    = 6;
 ops[2].parent_scaler_index = 2;
 ops[2].child1_clv_index    = 4;
 ops[2].child2_clv_index    = 5;
 ops[2].child1_matrix_index = 1;
 ops[2].child2_matrix_index = 0;
 ops[2].child1_scaler_index = 0;
 ops[2].child2_scaler_index = 1;

 /* use the operations array to compute 3 CLVs. Operations will be carried out
     starting from operation 0 to 2 */
 pll_update_partials(tree, ops, 3);

 /* print out the CLVs at tip and inner nodes*/
 /*
 printf ("Tip 0: ");
 pll_show_clv(tree,0,PLL_SCALE_BUFFER_NONE,6);
 printf ("Tip 1: ");
 pll_show_clv(tree,1,PLL_SCALE_BUFFER_NONE,6);
 printf ("Tip 2: ");
 pll_show_clv(tree,2,PLL_SCALE_BUFFER_NONE,6);
 printf ("Tip 3: ");
 pll_show_clv(tree,3,PLL_SCALE_BUFFER_NONE,6);
 printf ("CLV 4: ");
 pll_show_clv(tree,4,0,6);
 printf ("CLV 5: ");
 pll_show_clv(tree,5,1,6);
 printf ("CLV 6: ");
 pll_show_clv(tree,6,2,6);
 */

 /* compute the likelihood at the root of the rooted tree by specifying the CLV
     index of the root CLV and the index of the frequency vector to be used */
 double logl = pll_compute_root_loglikelihood(tree,6,2,0);

 printf("Log-L: %f\n", logl);



	pll_partition_destroy(tree);
	free(ops);
}
