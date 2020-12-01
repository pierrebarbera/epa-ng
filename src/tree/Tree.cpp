#include "tree/Tree.hpp"

#include <cstdio>
#include <iostream>
#include <numeric>
#include <stdexcept>

#include "core/pll/epa_pll_util.hpp"
#include "core/pll/optimize.hpp"
#include "io/file_io.hpp"
#include "seq/Sequence.hpp"
#include "util/logging.hpp"
#include "util/stringify.hpp"

Tree::Tree( std::string const& tree_file,
            MSA const& msa,
            raxml::Model& model,
            Options const& options,
            Memory_Footprint const& footprint )
    : ref_msa_( msa )
    , model_( model )
    , options_( options )
{
  try {
    tree_ = utree_ptr(
        build_tree_from_file( tree_file, nums_, mapper_, options.preserve_rooting ),
        utree_destroy );
  } catch( std::invalid_argument& e ) {
    auto modelstring = model_.num_states() == 4 ? "GTRGAMMAX" : "PROTGAMMAGTRX";
    std::cout << e.what() << " Please resolve the tree fully.\n( e.g. raxml -g "
              << tree_file << " -m " << modelstring << " -n <name> -s <alignment> -p 1234 )"
              << std::endl;
    throw std::runtime_error{ "Aborting" };
  }

  if( ref_msa_.size() != nums_.tip_nodes ) {
    LOG_WARN << "The reference MSA and tree have differing number of taxa! "
             << ref_msa_.size() << " vs. " << nums_.tip_nodes;
  }

  memsave_ = Memory_Config( options.memsave, footprint, tree_.get() );
  assert( memsave_ ? memsave_.structs : true );

  partition_
      = partition_ptr( make_partition( model_,
                                       tree_.get(),
                                       nums_,
                                       ref_msa_.num_sites(),
                                       options_,
                                       memsave_.structs ),
                       pll_partition_destroy );

  locks_ = Mutex_List( partition_->tips + partition_->clv_buffers );

  link_tree_msa( tree_.get(),
                 partition_.get(),
                 model_,
                 ref_msa_,
                 nums_.tip_nodes );

  set_unique_clv_indices( get_root( tree_.get() ), nums_.tip_nodes );

  branch_id_ = get_branch_ids( tree_.get() );

  LOG_DBG << model_;
  LOG_DBG << "Tree length: " << sum_branch_lengths( tree_.get() );

  pll_unode_t* vroot = nullptr;
  if( not memsave_ ) {
    precompute_clvs( tree_.get(), partition_.get(), nums_ );
  } else {
    // compute the CLVs toward the first vroot of the memsave-traversal
    vroot = memsave_.structs.traversal( 0 );
    partial_compute_clvs( tree_.get(),
                          nums_,
                          memsave_.structs.subtree_sizes(),
                          vroot,
                          partition_.get() );
  }

  auto logl = this->ref_tree_logl( vroot );

  if( logl == -std::numeric_limits< double >::infinity() ) {
    throw std::runtime_error{ "Tree Log-Likelihood -INF!" };
  }

  LOG_DBG << "Reference tree log-likelihood: "
          << std::to_string( logl );
}

/**
  Constructs the structures from binary file.
*/
Tree::Tree( std::string const& bin_file,
            raxml::Model& model,
            Options const& options,
            Memory_Footprint const& footprint )
    : model_( model )
    , options_( options )
    , binary_( bin_file )
{
  options_.load_binary_mode = true;
  partition_ = partition_ptr( binary_.load_partition(), pll_partition_destroy );
  nums_      = Tree_Numbers( partition_->tips );
  tree_  = utree_ptr( binary_.load_utree( partition_->tips ), utree_destroy );
  locks_ = Mutex_List( partition_->tips + partition_->clv_buffers );

  // TODO this needs a major facelift to be able to interoperate with all the
  // new shit
  branch_id_ = get_branch_ids( tree_.get() );

  memsave_ = Memory_Config( options.memsave, footprint, tree_.get() );

  if( memsave_ ) {
    throw std::runtime_error{"Temporarily disabled loading from binary into memsave partition"};
    // set the virtual root to where the traversal through the tree will later
    // (during the parallelized placement) start
    tree_.get()->vroot = memsave_.structs.traversal( 0 );

    // compute the CLVs toward that root
    partial_compute_clvs( tree_.get(),
                          nums_,
                          memsave_.structs.subtree_sizes(),
                          memsave_.structs.traversal( 0 ),
                          partition_.get() );
  }

  raxml::assign( model_, partition_.get() );
  LOG_DBG << model_;
  LOG_DBG << "Tree length: " << sum_branch_lengths( tree_.get() );
  LOG_DBG << "Reference tree log-likelihood: "
          << std::to_string( this->ref_tree_logl() );
}

/**
  Ensures that associated clv, tipchars, scalers are allocated and ready on
  return. If they are not currently in memory, fetches them from file if in
  Binary mode, or recomputes them in memsave mode.
*/
void Tree::ensure_clv_loaded( pll_unode_t const* const node )
{
  auto const clv_index = node->clv_index;

  // prevent race condition from concurrent access to this function
  Scoped_Mutex lock_by_clv_id( locks_[ clv_index ] );

  auto const scaler_index = node->scaler_index;
  bool const use_tipchars = partition_->attributes & PLL_ATTRIB_PATTERN_TIP;

  if( clv_index >= partition_->tips + partition_->clv_buffers ) {
    throw std::runtime_error{ "Node index out of bounds" };
  }

  if( use_tipchars and clv_index < partition_->tips ) {
    auto clv_ptr = partition_->tipchars[ clv_index ];
    // dynamically load from disk if not in memory
    if( options_.load_binary_mode
        and clv_ptr == nullptr ) {
      binary_.load_tipchars( partition_.get(), clv_index );
      clv_ptr = partition_->tipchars[ clv_index ];
    }
  } else {
    auto clv_ptr = pll_get_clv_reading(partition_.get(), clv_index );
    if (clv_ptr == nullptr) {
      if( options_.load_binary_mode ) {
        // dynamically load from disk if not in memory
        binary_.load_clv( partition_.get(), clv_index );
        clv_ptr = pll_get_clv_reading(partition_.get(), clv_index );
      } else if( this->memsave() ) {
        // kick off the partial traversal clv computation
        partial_compute_clvs( tree_.get(),
                              nums_,
                              memsave_.structs.subtree_sizes(),
                              const_cast< pll_unode_t* >( node ),
                              partition_.get() );
      } else {
        throw std::runtime_error{"Could not fetch CLV."};
      }
    }
  }

  // dynamically load the scaler_index if needed
  if( options_.load_binary_mode ) {
    if( scaler_index != PLL_SCALE_BUFFER_NONE
        and partition_->scale_buffer[ scaler_index ] == nullptr ) {
          binary_.load_scaler( partition_.get(), scaler_index );
    } else if( partition_->scale_buffer[ clv_index ] == nullptr ) {
      // the weird case where we are trying to load the scaler of a clv_index
      // (if they are indeed linked 1 to 1) for a scaler equaling
      // SCALE_BUFFER_NONE. This basically is only needed if we expect the
      // scaling situation to change, otherwise it's a pointless overhead.
      // since in epa-ng, the scaling cannot change and the binary holds a
      // partition that is "set in stone", we don't do anything here
    }
  }
}

double Tree::ref_tree_logl( pll_unode_t* const vroot )
{
  std::vector< unsigned int > param_indices( partition_->rate_cats, 0 );
  auto const root =  vroot ? vroot : get_root( tree_.get() );
  // ensure clvs are there
  this->ensure_clv_loaded( root );
  this->ensure_clv_loaded( root->back );

  auto logl = pll_compute_edge_loglikelihood( partition_.get(),
                                              root->clv_index,
                                              root->scaler_index,
                                              root->back->clv_index,
                                              root->back->scaler_index,
                                              root->pmatrix_index,
                                              &param_indices[ 0 ],
                                              nullptr );

  return logl;
}
