#include <iostream>
#include <unistd.h>
#include <string>
#include <algorithm>
#include <chrono>

#include "mpihead.hpp"
#include "logging.hpp"
#include "epa.hpp"

using namespace std;

Log lgr;

static void print_help()
{
  cout << "EPA - Evolutionary Placement Algorithm" << endl << endl;
  cout << "USAGE: epa [options]" << endl << endl;
  cout << "OPTIONS:" << endl;
  cout << "  -h \tDisplay this page" << endl;
  cout << "  -t \tPath to Tree file" << endl;
  cout << "  -s \tPath to reference MSA file. May also include the query sequences." << endl;
  cout << "  -q \tPath to separate query MSA file. If none is provided, epa will assume" << endl;
  cout << "     \tquery reads are in the reference MSA file (-s)" << endl;
  cout << "  -w \tPath to working directory" << endl;
  cout << "  -r \tRanged heuristic: only consider portion of query sequence sites not flanked by gaps during insertion" << endl;
  cout << "     \t  DEFAULT: OFF" << endl;
  cout << "  -g \tTwo-phase heuristic, determination of candidate edges using accumulative threshold" << endl;
  cout << "     \t" << endl;
  cout << "  -G \tTwo-phase heuristic, determination of candidate edges by specified percentage of total edges" << endl;
  cout << "     \t" << endl;
  cout << "  -O \toptimize reference tree and model parameters" << endl;
  cout << "  -l \tspecify minimum likelihood weight below which a placement is discarded" << endl;
  cout << "     \t  DEFAULT: 0.01" << endl;
  cout << "  -L \tspecify accumulated likelihood weight after which further placements are discarded" << endl;
  cout << "     \t  DEFAULT: OFF" << endl;
  cout << "  -m \tSpecify model of nucleotide substitution" <<  endl;
  cout << "     \tGTR \tGeneralized time reversible (DEFAULT)" << endl;
  cout << "     \tJC69\tJukes-Cantor Model" << endl;
  cout << "     \tK80 \tKimura 80 Model" << endl;
}

static void inv(string msg)
{
  int mpi_rank = 0;
  MPI_COMM_RANK(MPI_COMM_WORLD, &mpi_rank);
  if (mpi_rank == 0)
  {
    if (msg.size())
    lgr << msg << endl;
    print_help();
    lgr.flush();
  }
  MPI_FINALIZE();
  exit(EXIT_FAILURE);
}

static void ensure_dir_has_slash(string& dir)
{
  if (dir.length() > 0 && dir.back() != '/')
    dir += "/";
}

int main(int argc, char** argv)
{
  MPI_INIT(&argc, &argv);

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
  string tree_file("");
  string reference_file("");
  string binary_file("");

  int c;
  while((c =  getopt(argc, argv, "hOBb:t:s:q:l:L:w:g:G:r")) != EOF)
  {
    switch (c)
    {
      case 't':
        tree_file += optarg;
        break;
      case 's':
        reference_file += optarg;
        break;
      case 'q':
        query_file += optarg;
        break;
      case 'w':
        work_dir += optarg;
        break;
      case 'l':
        options.support_threshold = stod(optarg);
        if (options.support_threshold < 0.0)
          inv("Support threshold cutoff too small! Range: [0,1)");
        if (options.support_threshold >= 1.0)
          inv("Support threshold cutoff too large! Range: [0,1)");
        break;
      case 'L':
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
      case 'g':
        if (options.prescoring)
          inv("-g cannot be used simultaneously with -G!");
        if (optarg)
          options.prescoring_threshold = stod(optarg);
        if (options.prescoring_threshold < 0.0)
          inv("Prescoring threshold cutoff too small! Range: [0,1]");
        if (options.prescoring_threshold > 1.0)
          inv("Prescoring threshold cutoff too large! Range: [0,1]");
        options.prescoring = true;
        break;
      case 'G':
        if (options.prescoring)
          inv("-g cannot be used simultaneously with -G!");
        if (optarg)
          options.prescoring_threshold = stod(optarg);
        else
          options.prescoring_threshold = 0.1;
        if (options.prescoring_threshold < 0.0)
          inv("Prescoring threshold cutoff too small! Range: [0,1]");
        if (options.prescoring_threshold > 1.0)
          inv("Prescoring threshold cutoff too large! Range: [0,1]");
        options.prescoring = true;
        options.prescoring_by_percentage = true;
        break;
      case 'O':
        options.opt_branches = true;
        options.opt_model = true;
        break;
      case 'r':
        options.ranged = true;
        break;
      case 'B':
        options.dump_binary_mode = true;
        break;
      case 'b':
        options.load_binary_mode = true;
        binary_file += optarg;
        break;
      case ':':
        inv("Invalid option.");
        break;
      }
  }

  //================================================================
  //============    EPA    =========================================
  //================================================================

  ensure_dir_has_slash(work_dir);

  #ifdef __MPI
  int local_rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &local_rank);
  bool to_cout = (local_rank == 0);
  lgr = Log(work_dir + std::to_string(local_rank) + ".epa_info.log", to_cout);
  #else
  lgr = Log(work_dir + "epa_info.log");
  #endif

  MSA ref_msa;
  if (reference_file.size())
    ref_msa = build_MSA_from_file(reference_file);

  Model model(model_id);

  // build the Tree
  Tree tree;
  if (options.load_binary_mode)
  {
    lgr << "Loading from binary" << endl;
    tree = Tree(binary_file, options);
  }
  else
  {
    // build the full tree with all possible clv's
    tree = Tree(tree_file, ref_msa, model, options);
  }

  // build the query stream
  MSA_Stream queries;
  if (not options.dump_binary_mode)
  {
    if (query_file.size() != 0)
    {
      queries = MSA_Stream(query_file);

    }
    // attempt to split msa if it is intermingled with (supposed) query sequences
    else
    {
      throw runtime_error{"Combined MSA files not currently supported, please split them and specify using -s and -q."};
      // split_combined_msa(ref_msa, queries, tree);
    }
  }
  // dump to binary if specified
  if (options.dump_binary_mode)
  {
    lgr << "Writing to binary" << endl;
    string dump_file(work_dir + "epa_binary_file");
    dump_to_binary(tree, dump_file);
    MPI_FINALIZE();
  	return EXIT_SUCCESS;
  }

  // start the placement process and write to file
  auto start = chrono::high_resolution_clock::now();
  process(tree, queries, work_dir, options, invocation);
  auto end = chrono::high_resolution_clock::now();
  auto runtime = chrono::duration_cast<chrono::seconds>(end - start).count();

  lgr << "\nTime spent placing: " << runtime << "s" << endl;

  MPI_FINALIZE();
	return EXIT_SUCCESS;
}
