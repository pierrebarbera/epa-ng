#include <iostream>
#include <unistd.h>
#include <string>

#include "epa.hpp"

using namespace std;

static void print_help()
{
  cout << "EPA - Evolutionary PQuery Algorithm" << endl << endl;
  cout << "USAGE: epa [options] <reference_tree_file> <reference_MSA_file>" << endl << endl;
  cout << "OPTIONS:" << endl;
  cout << "  -h\t Display this page" << endl;
  cout << "  -q\t Path to separate query MSA file. If none is provided, epa will assume" << endl;
  cout << "    \t query reads are in the reference_MSA_file (second parameter)" << endl;
  cout << "  -b\t optimize branch lengths on insertion" << endl;
};

int main(int argc, char** argv)
{
  string invocation("");
  bool heuristic = true;
  for (int i = 0; i < argc; ++i)
  {
    invocation += argv[i];
    invocation += " ";
  }

  string query_file("");

  int c;
  while((c =  getopt(argc, argv, "hbq:")) != EOF)
  {
      switch (c)
      {
           case 'q':
               query_file += optarg;
               break;
           case 'h':
               print_help();
               exit(0);
               break;
           case 'b':
               heuristic = false;
               break;
           case ':':
               cerr << "Missing option." << endl;
               exit(EXIT_FAILURE);
               break;
      }
  }

  if (argc < 2)
  {
    cerr << "Insufficient parameters!" << endl;
    print_help();
    cout.flush();
    exit(EXIT_FAILURE);
  }
  // first two params are always the reference tree and msa file paths
  string tree_file(argv[optind]);
  string reference_file(argv[optind + 1]);

	epa(tree_file,
      reference_file,
      query_file,
      {0.25, 0.25, 0.25, 0.25},
      {1,1,1,1,1,1},
      1.0,
      heuristic,
      invocation);
	return 0;
}
