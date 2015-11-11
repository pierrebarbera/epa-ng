#ifndef EPA_FILE_IO_H_
#define EPA_FILE_IO_H_

#include <string>

class MSA;

/* reads in sequences from a file
  msa_file: string specifying the file path
*/
MSA* build_MSA_from_file(const std::string& msa_file);

#endif
