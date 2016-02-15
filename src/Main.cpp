#include <iostream>
#include <unistd.h>
#include <string>
#include <algorithm>

#include "epa.hpp"
#include "logging.hpp"

using namespace std;

static void print_help()
{
  cout << "EPA - Evolutionary Placement Algorithm" << endl << endl;
  cout << "USAGE: epa [options] <tree_file> <MSA_file>" << endl << endl;
  cout << "OPTIONS:" << endl;
  cout << "  -h \tDisplay this page" << endl;
  cout << "  -q \tPath to separate query MSA file. If none is provided, epa will assume" << endl;
  cout << "     \tquery reads are in the MSA_file (second parameter)" << endl;
  cout << "  -w \tPath to working directory" << endl;
  cout << "  -o \toptimize branch lengths on insertion" << endl;
  cout << "  -O \toptimize reference tree and model parameters" << endl;
  cout << "  -s \tspecify minimum likelihood weight below which a placement is discarded" << endl;
  cout << "     \t  DEFAULT: 0.01" << endl;
  cout << "  -S \tspecify accumulated likelihood weight after which further placements are discarded" << endl;
  cout << "     \t  DEFAULT: OFF" << endl;
  cout << "  -m \tSpecify model of nucleotide substitution" <<  endl;
  cout << "     \tGTR \tGeneralized time reversible`(DEFAULT)" << endl;
  cout << "     \tJC69 \tJukes-Cantor Model" << endl;
  cout << "     \tK80 \tKimura 80 Model" << endl;
}

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
  string work_dir("");

  int c;
  while((c =  getopt(argc, argv, "hoOq:s:S:w:")) != EOF)
  {
    switch (c)
    {
      case 'q':
        query_file += optarg;
        break;
      case 'w':
        work_dir += optarg;
        break;
      case 's':
        options.support_threshold = stod(optarg);
        if (options.support_threshold < 0.0)
          inv("Support threshold cutoff too small! Range: [0,1)");
        if (options.support_threshold >= 1.0)
          inv("Support threshold cutoff too large! Range: [0,1)");
        break;
      case 'S':
        options.support_threshold = stod(optarg);
        if (options.support_threshold <= 0.0)
          inv("Accumulated support threshold cutoff too small! Range: (0,1]");
        if (options.support_threshold > 1.0)
          inv("Accumulated support threshold cutoff too large! Range: (0,1]");
        options.acc_threshold = true;
        break;
      case 'h':
          inv("");
          break;
      case 'o':
          options.prescoring = true;
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
    work_dir,
    model,
    options,
    invocation);

	return 0;
}
