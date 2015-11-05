#include "MSA.h"

using namespace std;

MSA::MSA(const string& msa_file)
{
  build_from_file(msa_file);
}

MSA::~MSA()
{
  // vectors are self deleting
}

/* reads in sequences from a file
  msa_file: string specifying the file path
*/
void MSA::build_from_file(const string& msa_file)
{

  /* open the file */
  auto file = pll_fasta_open(msa_file.c_str(), pll_map_fasta);
  if (!file)
    throw runtime_error{string("Cannot open file") + msa_file};

  char * sequence = NULL;
  char * header = NULL;
  long sequence_length;
  long header_length;
  long sequence_number;

  /* read sequences and make sure they are all of the same length */
  int num_sites = 0;
  
  while (pll_fasta_getnext(file, &header, &header_length, &sequence, &sequence_length, &sequence_number))
  {

    if (num_sites && num_sites != sequence_length)
      throw runtime_error{"MSA file does not contain equal size sequences"};

    if (!num_sites) num_sites = sequence_length;

    headers.push_back(header);
    sequences.push_back(sequence);
  }

  if (pll_errno != PLL_ERROR_FILE_EOF)
    throw runtime_error{string("Error while reading file:") +  msa_file};

  pll_fasta_close(file);

  if (num_sites == -1)
    throw runtime_error{"Unable to read MSA file"};

  MSA::num_sites = num_sites;


}

tuple<string, string> MSA::get(int i)
{
  return make_tuple(headers[i], sequences[i]);
}
