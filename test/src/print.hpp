#pragma once

#include "core/pll/pllhead.hpp"
#include "io/file_io.hpp"
#include "seq/MSA.hpp"
#include "seq/MSA_Info.hpp"
#include "tree/Tiny_Tree.hpp"
#include "tree/Tree.hpp"
#include "util/Options.hpp"

inline void print_tipchars( pll_partition_t const& part)
{
  if( part.tipchars ) {
    for( size_t tip_id = 0; tip_id < part.tips; tip_id++ ) {
      std::cout << "tip " << tip_id << "\n";
      for( size_t site = 0; site < part.sites; site++ ) {
        if( part.tipchars[ tip_id ] ) {
          std::cout << (char)part.tipchars[ tip_id ][ site ];
        }
      }
      std::cout << "\n";
    }
  } else {
    std::cout << "No tipchars array allocd\n";
  }
}