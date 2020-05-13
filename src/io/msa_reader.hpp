#pragma once

#include <memory>

#include "io/Binary_Fasta.hpp"
#include "io/file_io.hpp"
#include "io/msa_reader_interface.hpp"
#include "seq/MSA_Info.hpp"
#include "seq/MSA_Stream.hpp"
#include "util/Options.hpp"
#include "util/logging.hpp"
#include "util/stringify.hpp"

inline auto make_msa_reader( std::string const& file_name,
                             MSA_Info const& info,
                             bool const premasking = true,
                             bool const split      = false )
{
  std::unique_ptr< msa_reader > result( nullptr );

  try {
    result = std::make_unique< Binary_Fasta_Reader >( file_name, info, premasking, split );
  } catch( std::exception const& e ) {
    LOG_DBG << "Failed to parse input as binary fasta (bfast), trying `fasta` instead.";
    result = std::make_unique< MSA_Stream >( file_name, info, premasking, split );
  }

  return result;
}
