#include "genesis/genesis.hpp"

#include "util/logging.hpp"

static genesis::sequence::SequenceSet read_any_seqfile( std::string const& file )
{
  genesis::sequence::SequenceSet out_set;

  auto reader = genesis::sequence::PhylipReader();
  try {
    reader.read( genesis::utils::from_file( file ), out_set );
  } catch( std::exception& e ) {
    out_set.clear();
    reader.mode( genesis::sequence::PhylipReader::Mode::kInterleaved );
    try {
      reader.read( genesis::utils::from_file( file ), out_set );
    } catch( std::exception& e ) {
      out_set.clear();
      try {
        genesis::sequence::FastaReader().read( genesis::utils::from_file( file ), out_set );
      } catch( std::exception& e ) {
        throw std::invalid_argument{
          "Cannot parse sequence file(s): Invalid file format? (only phylip and fasta allowed)"
        };
      }
    }
  }

  return out_set;
}

static void write_subset( genesis::sequence::SequenceSet const& set,
                          std::unordered_set< std::string > subset_labels,
                          std::string const& output_file )
{
  genesis::sequence::SequenceSet ref_set;

  // for every sequence in the qry, until all ref are found
  for( auto const& seq : set ) {
    auto it = subset_labels.find( seq.label() );

    if( it != subset_labels.end() ) {
      // add the sequence to the reference
      ref_set.add( seq );
      // shrink the set of labels
      subset_labels.erase( it );
    }

    // contine while there are still labels
    if( subset_labels.empty() ) {
      break;
    }
  }

  // ensure all references were found
  if( not subset_labels.empty() ) {
    throw std::invalid_argument{
      "Could not find all references in the query file to produce an appropriate reference.fasta!"
    };
  }

  // output the ref file
  genesis::sequence::FastaWriter().to_file( ref_set, output_file );
}

void split( std::string const& ref_msa,
            std::vector< std::string > query_files,
            std::string const& outdir = "" )
{
  auto outfile = outdir + "query.fasta";

  auto ref_set           = read_any_seqfile( ref_msa );
  auto ref_labels        = labels( ref_set );
  size_t const ref_width = ref_set.at( 0 ).length();

  std::ofstream outstream;
  genesis::utils::file_output_stream( outfile, outstream );

  size_t qry_width     = 0;
  bool first_iteration = true;
  auto writer          = genesis::sequence::FastaWriter();
  for( auto const& qry_file : query_files ) {

    auto qry_set               = read_any_seqfile( qry_file );
    size_t const cur_qry_width = qry_set.at( 0 ).length();

    // check for consistent query width across files
    if( qry_width and cur_qry_width != qry_width ) {
      throw std::invalid_argument{ std::string( "Query file '" ) + qry_file + "' appears to have an alignment width"
                                   + " that differs from previous query alignments. Aborting!" };
    }

    qry_width = cur_qry_width;

    // peek into query file to check if the alignment width agrees with the reference
    if( cur_qry_width != ref_width ) {
      LOG_WARN << "The query alignment file '" << qry_file
               << "' appears to have an alignment width that differs from the reference ("
               << cur_qry_width << " vs. " << ref_width << ").\n"
               << "This is likely due to the alignment tool stripping gap-only columns, or adding columns to the reference."
               << " Please consider using the produced 'reference.fasta' during placement!";
    }

    if( first_iteration ) {
      // write the appropriate ref file
      write_subset( qry_set, ref_labels, outdir + "reference.fasta" );
      first_iteration = false;
    }

    genesis::sequence::filter_by_label_list( qry_set, ref_labels );

    writer.to_stream( qry_set, outstream );
  }
}
