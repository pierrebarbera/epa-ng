#include "genesis/genesis.hpp"

#include "util/logging.hpp"

static genesis::sequence::SequenceSet read_any_seqfile(std::string const& file)
{
    genesis::sequence::SequenceSet out_set;

    auto reader = genesis::sequence::PhylipReader();
    try {
        reader.read( genesis::utils::from_file( file ), out_set );
    } catch(std::exception& e) {
        out_set.clear();
        reader.mode( genesis::sequence::PhylipReader::Mode::kInterleaved );
        try {
            reader.read( genesis::utils::from_file( file ), out_set );
        } catch(std::exception& e) {
            out_set.clear();
            try {
                genesis::sequence::FastaReader().read( genesis::utils::from_file( file ), out_set );
            } catch(std::exception& e) {
                throw std::invalid_argument{"Cannot parse sequence file(s): Invalid file format? (only phylip and fasta allowed)"};
            }
        }
    }

    return out_set;
}

void split(std::string ref_msa, std::vector<std::string> query_files, std::string outdir="")
{
    auto outfile = outdir + "query.fasta";

    auto ref_set = read_any_seqfile(ref_msa);
    auto ref_labels = labels( ref_set );


    std::ofstream outstream;
    genesis::utils::file_output_stream(outfile, outstream);

    auto writer = genesis::sequence::FastaWriter();
    for (auto const& qry_file : query_files) {

        auto qry_set = read_any_seqfile(qry_file);
        genesis::sequence::filter_by_label_list( qry_set, ref_labels );

        writer.to_stream( qry_set, outstream );
    }
}
