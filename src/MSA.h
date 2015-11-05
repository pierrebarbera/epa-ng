#include <vector>
#include <string>
#include <assert.h>
#include <stdexcept>

#ifndef __PLL__
#define __PLL__
#include "pll.h"
#endif


class MSA
{
public:
  MSA(const std::string& msa_file);
  ~MSA();
  void build_from_file(const std::string& msa_file);

  int num_sites;

private:

  std::vector<std::string>* headers;
  std::vector<std::string>* sequences;
};