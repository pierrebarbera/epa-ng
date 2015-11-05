#include "MSA.h"

MSA::MSA(const std::string& msa_file)
{
  this->headers = new std::vector<std::string>;
  this->sequences = new std::vector<std::string>;
  build_from_file(msa_file);
}

MSA::~MSA()
{
  delete this->headers;
  delete this->sequences;
}

/* reads in sequences from a file
  msa_file: string specifying the file path
*/
void MSA::build_from_file(const std::string& msa_file)
{

  /* open the file */
  auto file = pll_fasta_open(msa_file.c_str(), pll_map_fasta);
  if (!file)
    throw std::runtime_error{std::string("Cannot open file") + msa_file};

  char * sequence = NULL;
  char * header = NULL;
  long sequence_length;
  long header_length;
  long sequence_number;

  /* read sequences and make sure they are all of the same length */
  int num_sites = 0;
  
  while (pll_fasta_getnext(file, &header, &header_length, &sequence, &sequence_length, &sequence_number))
  {

    if (num_sites != sequence_length)
      throw std::runtime_error{"MSA file does not contain equal size sequences"};

    if (!num_sites) num_sites = sequence_length;

    this->headers->push_back(header);
    this->sequences->push_back(sequence);
  }

  if (pll_errno != PLL_ERROR_FILE_EOF)
    throw std::runtime_error{std::string("Error while reading file:") +  msa_file};

  pll_fasta_close(file);

  if (num_sites == -1)
    throw std::runtime_error{"Unable to read MSA file"};

  this->num_sites = num_sites;


}
