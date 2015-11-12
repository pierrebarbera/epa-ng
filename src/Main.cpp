#include <iostream>
#include <string>

#include "epa.hpp"

using namespace std;

void plltest();

int main(int argc, char** argv)
{
  string tree_file = "";
  string reference_file = "";
  string query_file = "";

  if(argc != 0)
  {
    tree_file = string(argv[1]);
    reference_file = string(argv[2]);
    query_file = string(argv[3]);
  } else {
    exit(EXIT_FAILURE);
  }

	//plltest();

  //auto tree_file = new string("abc");
	//auto msa_file = new string("def");
	epa(tree_file,
      reference_file,
      query_file,
      {0.25, 0.25, 0.25, 0.25},
      {1,1,1,1,1,1});
	// status callback?
	printf("gday\n");
 return 0;
}
