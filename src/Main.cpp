#include <iostream>
#include <unistd.h>
#include <string>

#include "epa.hpp"

using namespace std;

static void print_help()
{
  cout << "EPA - Evolutionary Placement Algorithm" << endl << endl;
  cout << "USAGE: epa [options] <reference_tree_file> <reference_MSA_file>" << endl << endl;
  cout << "OPTIONS:" << endl;
  cout << "  -h \tDisplay this page" << endl;
  cout << "  -q \tPath to separate query MSA file. If none is provided, epa will assume" << endl;
  cout << "     \tquery reads are in the reference_MSA_file (second parameter)" << endl;
  cout << "  -o \toptimize branch lengths on insertion" << endl;
  cout << "  -O \toptimize reference tree and model parameters" << endl;
  cout << "  -m \tSpecify model of nucleotide substitution" << endl << endl;
  cout << "     \tGTR \tGeneralized time reversible" << endl;
  cout << "     \tJC69 \tJukes-Cantor Model" << endl;
  cout << "     \tK80 \tKimura 80 Model" << endl;
};

static void inv(string msg)
{
  if (msg.size())
    cerr << msg << endl;
  print_help();
  cout.flush();
  exit(EXIT_FAILURE);
}

int main(int argc, char** argv)
{
  string invocation("");
  string model_id("GTR");
  Options options;
  for (int i = 0; i < argc; ++i)
  {
    invocation += argv[i];
    invocation += " ";
  }

  string query_file("");

  int c;
  while((c =  getopt(argc, argv, "hoOq:s:")) != EOF)
  {
    switch (c)
    {
      case 'q':
        query_file += optarg;
        break;
      case 's':
        options.support_threshold = stod(optarg);
        if (options.support_threshold <= 0.0)
          inv("Support threshold cutoff too small!");
        break;
      case 'h':
          inv("");
          break;
      case 'o':
          options.opt_insertion_branches = true;
          break;
      case 'O':
          options.opt_branches = true;
          options.opt_model = true;
          break;
      case ':':
          inv("Missing option.");
          break;
      }
  }

  if (argc < 2)
    inv("Insufficient parameters!");

  // first two params are always the reference tree and msa file paths
  string tree_file(argv[optind]);
  string reference_file(argv[optind + 1]);

  Model model(model_id);

	epa(tree_file,
      reference_file,
      query_file,
      model,
      options,
      invocation);
	return 0;
}
