#include "seq/MSA_Info.hpp"

#include "io/Binary_Fasta.hpp"

MSA_Info make_msa_info( std::string const& file_path )
{
  MSA_Info info;
  try {
    info = Binary_Fasta::get_info( file_path );
  } catch( std::exception const& ) {
    info = MSA_Info( file_path );
  }
  return info;
}