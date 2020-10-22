#pragma once

#include <array>
#include <cassert>
#include <limits>
#include <map>
#include <mutex>
#include <utility>
#include <vector>

#include "seq/Sequence.hpp"
#include "tree/Tiny_Tree.hpp"
#include "util/Matrix.hpp"
#include "util/Range.hpp"
#include "util/maps.hpp"

constexpr size_t INVALID = std::numeric_limits< size_t >::max();

class Lookup_Store {
  /**
   * NOTE TO FUTURE DEVS:
   * This class has gotten a bit convoluted, so where is a brief overview of the
   * various maps and tables:
   *
   * store_: vector of matrices, one per branch in the ref tree
   *
   * <matrix in store>-> lookup_matrix: stores one CLV per character suitable
   * for the model (ACGTVH- etc.)
   *
   * char_map_: set of chars for which a lookup_matrix is done (see
   * util/maps.hpp)
   *
   * char_to_posish: maps ascii char to a column in a lookup_matrix. This also
   * normalizes the input! meaning: map upper and lowercase to the same CLV
   * site, different variants of GAP (-?Xx etc.) and ANY (N), U into T
   * (RNA support) and defines invalid chars
   */
  public:
  using lookup_type = Matrix< double >;

  Lookup_Store( size_t const num_branches, size_t const num_states )
      : branch_( num_branches )
      , store_( num_branches )
      , char_map_size_( ( num_states == 4 ) ? NT_MAP_SIZE : AA_MAP_SIZE )
      , char_map_( ( num_states == 4 ) ? NT_MAP : AA_MAP )
  {
    bool const dna = ( num_states == 4 );

    for( size_t i = 0; i < 128; ++i ) {
      char_to_posish_[ i ] = INVALID;
    }

    // build reverse map from char to ID in the charmap
    for( size_t i = 0; i < char_map_size_; ++i ) {
      char_to_posish_[ char_map_[ i ] ]                 = i;
      char_to_posish_[ std::tolower( char_map_[ i ] ) ] = i;
    }

    if( dna ) {
      // allow for RNA
      char_to_posish_[ 'U' ] = char_to_posish_[ 'T' ];
      char_to_posish_[ 'u' ] = char_to_posish_[ 'T' ];
    }

    // gap/any chars
    if( dna ) {
      char_to_posish_[ 'X' ] = char_to_posish_[ '-' ];
      char_to_posish_[ 'x' ] = char_to_posish_[ '-' ];
      char_to_posish_[ 'O' ] = char_to_posish_[ '-' ];
      char_to_posish_[ 'o' ] = char_to_posish_[ '-' ];
      char_to_posish_[ '.' ] = char_to_posish_[ '-' ];
    } else {
      char_to_posish_[ 'X' ] = char_to_posish_[ 'N' ];
      char_to_posish_[ 'x' ] = char_to_posish_[ 'N' ];
    }
    char_to_posish_[ '?' ] = char_to_posish_[ '-' ];
  }

  Lookup_Store()  = delete;
  ~Lookup_Store() = default;

  void init_branch( size_t const branch_id,
                    std::vector< std::vector< double > > precomps )
  {
    store_[ branch_id ]
        = Matrix< double >( precomps[ 0 ].size(), char_map_size_ );

    for( size_t ch = 0; ch < precomps.size(); ++ch ) {
      for( size_t site = 0; site < precomps[ ch ].size(); ++site ) {
        store_[ branch_id ]( site, ch ) = precomps[ ch ][ site ];
      }
    }
  }

  void init_branch( Tiny_Tree const& tiny_tree )
  {
    auto const size = char_map_size();

    // precompute all possible site likelihoods
    std::vector< std::vector< double > > precomputed_sites( size );
    for( size_t i = 0; i < size; ++i ) {
      tiny_tree.get_persite_logl( char_map( i ), precomputed_sites[ i ] );
    }
    init_branch( tiny_tree.branch_id(), precomputed_sites );
  }

  std::mutex& get_mutex( size_t const branch_id )
  {
    return branch_[ branch_id ];
  }

  bool has_branch( size_t const branch_id )
  {
    return store_[ branch_id ].size() != 0;
  }

  lookup_type& operator[]( size_t const branch_id )
  {
    return store_[ branch_id ];
  }

  unsigned char char_map( size_t const i )
  {
    if( i >= char_map_size_ ) {
      throw std::runtime_error{ std::string(
                                    "char_map access out of bounds! i =" )
                                + std::to_string( i ) };
    }

    return char_map_[ i ];
  }

  size_t char_map_size() { return char_map_size_; }

  size_t char_position( unsigned char c ) const
  {
    auto pos = char_to_posish_[ c ];

    if( pos == INVALID ) {
      throw std::runtime_error{ std::string( "char is invalid! char = " )
                                + std::to_string( c ) };
    }

    return pos;
  }

  /**
   * Sum up the per-site log-likelihoods for a given branch, based on the given
   * sequence. This LH is dependant on the values that the Tiny_Tree had when
   * it was used to create the lookup. The only real use case is preplacement,
   * meaning placed in the center of the insertion branch, with the pednant
   * length at the default value.
   *
   * @param  branch_id  branch for which the liklihood should be calculated
   * @param  seq        Sequence for which likelihood is to be calculated
   * @param  premasking flag indicating if flanking gaps should be masked out
   * @return            log-likelihood
   */
  double sum_precomputed_sitelk( size_t const branch_id,
                                 std::string const& seq,
                                 bool const premasking ) const
  {
    Range range( 0, seq.size() );
    if( premasking ) {
      range = get_valid_range( seq );
      if( not range ) {
        throw std::runtime_error{
          "A sequence does not appear to have any non-gap sites!"
        };
      }
    }

    assert( seq.length() == store_[ branch_id ].rows() );

    double sum                = 0;
    auto const& lookup_matrix = store_[ branch_id ];
    auto const& lookup        = lookup_matrix.get_array();

    // unrolled loop
    size_t site      = range.begin;
    size_t const end = range.begin + range.span;

    size_t const stride = 4;
    for( ; site + stride - 1u < end; site += stride ) {
      double sum_one = lookup[ lookup_matrix.coord(
                           site, char_to_posish_[ seq[ site ] ] ) ]
          + lookup[ lookup_matrix.coord(
                site + 1u, char_to_posish_[ seq[ site + 1u ] ] ) ];

      double sum_two = lookup[ lookup_matrix.coord(
                           site + 2u, char_to_posish_[ seq[ site + 2u ] ] ) ]
          + lookup[ lookup_matrix.coord(
                site + 3u, char_to_posish_[ seq[ site + 3u ] ] ) ];

      sum_one += sum_two;

      sum += sum_one;
    }

    // rest of the horizontal add
    while( site < end ) {
      sum += lookup[ lookup_matrix.coord( site,
                                          char_to_posish_[ seq[ site ] ] ) ];
      ++site;
    }
    return sum;
  }

  private:
  std::vector< std::mutex > branch_;
  std::vector< lookup_type > store_;
  size_t const char_map_size_;
  unsigned char const* char_map_;
  std::array< size_t, 128 > char_to_posish_;
};

#ifdef __OMP
#include <omp.h>
#endif

#include "core/pll/pllhead.hpp"
#include "core/BranchBuffer.hpp"
#include "tree/Tree.hpp"

/**
 * Wraps the lookup store such that this object can be queried for a
 * lookup-based
 * placement (used for preplacement) for a given QS on a given branch.
 *
 * Creation of the object should take care of creating all necessary loopups,
 * do it nicely parallel, and ideally optionally using the memsaver.
 */
class LookupPlacement {
  public:
  LookupPlacement( Tree& ref_tree,
                   std::vector< pll_unode_t* > const& branches,
                   Options const& options )
      : lookup_( ref_tree.nums().branches, ref_tree.partition()->states )
      , pendant_length_( ref_tree.nums().branches, -1.0 )
      , distal_length_( ref_tree.nums().branches, -1.0 )
  {
#ifdef __OMP
    omp_set_num_threads( options.num_threads ? options.num_threads
                                             : omp_get_max_threads() );
#endif
    auto nums = ref_tree.nums();
    bool const use_memsave
        = ( ref_tree.partition()->attributes & PLL_ATTRIB_LIMIT_MEMORY );

    // create and hold the lookup table for the entirety of the reference tree

    // if the partition is not in memsave mode, calculate the lookups normally,
    // from the fully precomputed partition data
    if( not use_memsave ) {
#ifdef __OMP
#pragma omp parallel for schedule( dynamic )
#endif
      for( size_t branch_id = 0; branch_id < nums.branches; ++branch_id ) {
        Tiny_Tree cur_branch( branches[ branch_id ], branch_id, ref_tree );

        pendant_length_[ branch_id ] = cur_branch.pendant_length();
        distal_length_[ branch_id ]  = cur_branch.distal_length();

        lookup_.init_branch( cur_branch );
      }

    } else {
      // otherwise do it in a regulated way ensuring limited memory use
      auto const block_size = options.memory_config.concurrent_branches;
      BranchBuffer branchbuf( &ref_tree, block_size );
      BranchBuffer::container_type branch_chunk;

      while( branchbuf.get_next( branch_chunk ) ) {
        // parallelize over branches: each thread places all queries on its
        // designated branch
#pragma omp parallel for schedule( dynamic )
        for( size_t i = 0; i < branch_chunk.size(); ++i ) {
          auto& cur_branch     = branch_chunk[ i ];
          auto const branch_id = cur_branch.branch_id();

          pendant_length_[ branch_id ] = cur_branch.pendant_length();
          distal_length_[ branch_id ]  = cur_branch.distal_length();

          lookup_.init_branch( cur_branch );
        }
      }
    }
  }

  // deleted copy assignment since this class hold hella memory

  LookupPlacement()  = delete;
  ~LookupPlacement() = default;

  LookupPlacement( LookupPlacement const& other ) = delete;
  LookupPlacement( LookupPlacement&& other )      = default;

  LookupPlacement& operator=( LookupPlacement const& other ) = delete;
  LookupPlacement& operator=( LookupPlacement&& other ) = default;

  Placement place( size_t const branch_id,
                   Sequence const& seq,
                   bool const premasking ) const
  {
    auto const logl = lookup_.sum_precomputed_sitelk(
        branch_id, seq.sequence(), premasking );

    return Placement( branch_id,
                      logl,
                      pendant_length_[ branch_id ],
                      distal_length_[ branch_id ] );
  }

  size_t num_branches() const { return pendant_length_.size(); }

  private:
  Lookup_Store lookup_;
  // arrays holding pendant and distal lengths for a given branch_id
  std::vector< double > pendant_length_;
  std::vector< double > distal_length_;
};
