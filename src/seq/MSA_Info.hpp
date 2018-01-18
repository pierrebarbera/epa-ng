#include "genesis/sequence/functions/functions.hpp"
#include "genesis/sequence/formats/fasta_input_iterator.hpp"

#include <string>

/**
 * Class encompassing info about a MSA File.
 *
 * Mainly used for storing masking and size info.
 */
class MSA_Info
{
  using Bitvector = genesis::utils::Bitvector;

public:

  /**
   * Pass through the given file and generate info
   */
  MSA_Info( const std::string& file_path,
            std::function<void(const genesis::sequence::Sequence&)> fn = nullptr )
    : path_(file_path)
  {
    // detect number of sequences in fasta file and generate mask
    genesis::utils::InputStream instr( std::make_unique< genesis::utils::FileInputSource >( file_path ));
    auto it = genesis::sequence::FastaInputIterator( instr );

    // set some initial stuff
    if (it) {
      sites_ = it->length();
      gap_mask_ = Bitvector(sites_, true);
    }

    std::vector<size_t> entry_sizes;

    while ( it ) {
      const auto& seq = *it;

      // call the supplied lambda, if valid
      if (fn) {
        fn(seq);
      }

      // ensure sequence have consistent length
      if (sites_ and (sites_ != seq.length())) {
        throw std::runtime_error{path_
          + " does not contain equal size sequences! First offending sequence: "
          + seq.label()};
      }

      // get the mask of the current sequence
      auto cur_mask = genesis::sequence::gap_sites(seq);
      // adjust global mask accordingly
      gap_mask_ &= cur_mask;

      ++it;
    }
  }
  MSA_Info() = default;
  ~MSA_Info() = default;

  // access
  const std::string& path() const {return path_;}
  size_t sites() const {return sites_;}
  const Bitvector& gap_mask() const {return gap_mask_;}
  size_t gap_count() const {return gap_mask_.count();}

private:
  std::string path_ = "";
  size_t sites_ = 0;
  Bitvector gap_mask_;
  
};

inline std::ostream& operator << (std::ostream& out, MSA_Info const& rhs)
{
  out << "Path: " << rhs.path();
  out << "\nSites: " << rhs.sites();
  out << "\nGaps: " << rhs.gap_count();
  out << "\nFraction of gaps: " << rhs.gap_count() / static_cast<double>(rhs.sites());
  out << "\n";

  return out;
}
