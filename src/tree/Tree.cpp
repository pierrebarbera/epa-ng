#include "tree/Tree.hpp"

#include <cstdio>
#include <iostream>
#include <numeric>
#include <stdexcept>

#include "core/pll/epa_pll_util.hpp"
#include "core/pll/optimize.hpp"
#include "io/file_io.hpp"
#include "seq/Sequence.hpp"
#include "set_manipulators.hpp"
#include "util/logging.hpp"
#include "util/stringify.hpp"

Tree::Tree( std::string const& tree_file,
            MSA const& msa,
            raxml::Model& model,
            Options const& options )
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

  if( options_.memsave ) {
    memsave_ = Memsaver( tree_.get() );
  }

  partition_ = partition_ptr( make_partition( model_,
                                              nums_,
                                              ref_msa_.num_sites(),
                                              options_,
                                              subtree_sizes_ ),
                              pll_partition_destroy );

  locks_ = Mutex_List( partition_->tips + partition_->clv_buffers );

  link_tree_msa( tree_.get(),
                 partition_.get(),
                 model_,
                 ref_msa_,
                 nums_.tip_nodes );

  set_unique_clv_indices( get_root( tree_.get() ), nums_.tip_nodes );

  branch_id_ = get_branch_ids( tree_ );

  LOG_DBG << model_;
  LOG_DBG << "Tree length: " << sum_branch_lengths( tree_.get() );

  if( not memsave_ ) {
    precompute_clvs( tree_.get(), partition_.get(), nums_ );
  } else {
    // set the virtual root to where the traversal through the tree will later
    // (during the parallelized placement) start
    tree_.get()->vroot = memsave_.traversal( 0 );

    // compute the CLVs toward that root
    partial_compute_clvs( tree_.get(),
                          nums_,
                          subtree_sizes_.get(),
                          traversal[ 0 ],
                          partition_.get() );
  }

  auto logl = this->ref_tree_logl();

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
            Options const& options )
    : model_( model )
    , options_( options )
    , binary_( bin_file )
{
  options_.load_binary_mode = true;
  partition_                = partition_ptr( binary_.load_partition(), pll_partition_destroy );
  nums_                     = Tree_Numbers( partition_->tips );
  tree_                     = utree_ptr( binary_.load_utree( partition_->tips ), utree_destroy );
  locks_                    = Mutex_List( partition_->tips + partition_->clv_buffers );

  // TODO this needs a major facelift to be able to interoperate with all the new shit
  branch_id_ = get_branch_ids( tree_ );
  if( options_.memsave ) {
    memsave_ = Memsaver( tree_.get() );
    // set the virtual root to where the traversal through the tree will later
    // (during the parallelized placement) start
    tree_.get()->vroot = memsave_.traversal( 0 );

    // compute the CLVs toward that root
    partial_compute_clvs( tree_.get(),
                          nums_,
                          subtree_sizes_.get(),
                          traversal[ 0 ],
                          partition_.get() );
  }

  raxml::assign( model_, partition_.get() );
  LOG_DBG << model_;
  LOG_DBG << "Tree length: " << sum_branch_lengths( tree_.get() );
  LOG_DBG << "Reference tree log-likelihood: "
          << std::to_string( this->ref_tree_logl() );
}

/**
  Returns a pointer either to the CLV or tipchar buffer, depending on the index.
  If they are not currently in memory, fetches them from file.
  Ensures that associated scalers are allocated and ready on return.
*/
void* Tree::get_clv( pll_unode_t const* const node )
{
  auto const i = node->clv_index;

  // prevent race condition from concurrent access to this function
  Scoped_Mutex lock_by_clv_id( locks_[ i ] );

  auto const scaler       = node->scaler_index;
  bool const use_tipchars = partition_->attributes & PLL_ATTRIB_PATTERN_TIP;

  if( i >= partition_->tips + partition_->clv_buffers ) {
    throw std::runtime_error{ "Node index out of bounds" };
  }

  void* clv_ptr = nullptr;
  if( use_tipchars and i < partition_->tips ) {
    clv_ptr = partition_->tipchars[ i ];
    // dynamically load from disk if not in memory
    if( options_.load_binary_mode
        and clv_ptr == nullptr ) {
      binary_.load_tipchars( partition_.get(), i );
      clv_ptr = partition_->tipchars[ i ];
    }
  } else {
    clv_ptr = partition_->clv[ i ];
    // dynamically load from disk if not in memory
    if( options_.load_binary_mode
        and clv_ptr == nullptr ) {
      binary_.load_clv( partition_.get(), i );
      clv_ptr = partition_->clv[ i ];
    }
  }

  // dynamically load the scaler if needed
  if( options_.load_binary_mode
      and scaler != PLL_SCALE_BUFFER_NONE
      and partition_->scale_buffer[ scaler ] == nullptr ) {
    binary_.load_scaler( partition_.get(), scaler );
  }

  return clv_ptr;
}

double Tree::ref_tree_logl()
{
  std::vector< unsigned int > param_indices( partition_->rate_cats, 0 );
  auto const root = get_root( tree_.get() );
  // ensure clvs are there
  this->get_clv( root );
  this->get_clv( root->back );

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
