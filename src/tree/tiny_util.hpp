#pragma once

#include "core/pll/pllhead.hpp"
#include "tree/Tree.hpp"

void tiny_partition_destroy( pll_partition_t* partition );
pll_utree_t* make_tiny_tree_structure( pll_unode_t const* old_proximal,
                                       pll_unode_t const* old_distal,
                                       bool const tip_tip_case );
pll_partition_t* make_tiny_partition( Tree& reference_tree,
                                      pll_utree_t const* tree,
                                      pll_unode_t const* old_proximal,
                                      pll_unode_t const* old_distal,
                                      bool const tip_tip_case );
