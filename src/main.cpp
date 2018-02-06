#include <iostream>
#include <string>
#include <algorithm>
#include <chrono>

#include <cxxopts.hpp>

#include "net/mpihead.hpp"
#include "util/logging.hpp"
#include "util/Options.hpp"
#include "util/stringify.hpp"
#include "util/parse_model.hpp"
#include "io/Binary_Fasta.hpp"
#include "io/Binary.hpp"
#include "io/file_io.hpp"
#include "io/msa_reader.hpp"
#include "tree/Tree.hpp"
#include "core/raxml/Model.hpp"
#include "core/place.hpp"
#include "seq/MSA.hpp"
#include "seq/MSA_Info.hpp"

static void ensure_dir_has_slash(std::string& dir)
{
  if (dir.length() > 0 && dir.back() != '/') {
    dir += "/";
  }
}

inline bool is_file (const std::string& name) {
    std::ifstream f(name.c_str());
    return f.good();
}

void exit_epa(int ret=EXIT_SUCCESS)
{
  MPI_FINALIZE();
  std::exit(ret);
}

int main(int argc, char** argv)
{
  auto start_all = std::chrono::high_resolution_clock::now();
  genesis::utils::Logging::log_to_stdout();

#ifdef __MPI
  MPI_INIT(&argc, &argv);
  int local_rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &local_rank);
  if (local_rank != 0) {
    genesis::utils::Logging::log_to_stdout(false);  
  }
#endif
  genesis::utils::Logging::max_level(genesis::utils::Logging::kInfo);

  std::string invocation("");
  std::string model_desc("GTR+G");
  Options options;

  for (int i = 0; i < argc; ++i) {
    invocation += argv[i];
    invocation += " ";
  }

  std::string query_file("");
  std::string work_dir(".");
  std::string tree_file("");
  std::string reference_file("");
  std::string binary_file("");

  std::string banner;

  raxml::Model model;

  const bool empty = argc == 1;
  bool pipeline = false;

  try
  {
  cxxopts::Options cli(argv[0], "Massively-Parallel Evolutionary Placement Algorithm");

  cli.add_options()
    ("help", "Display help.")
    ("v,version", "Display version.")
    ("verbose", "Display debug information.")
    ;
  cli.add_options("Input")
    ("t,tree", "Path to Reference Tree file.", cxxopts::value<std::string>())
    ("s,ref-msa", "Path to Reference MSA file.", cxxopts::value<std::string>())
    ("q,query", "Path to Query MSA file.", cxxopts::value<std::string>())
    ("b,binary", "Path to Binary file.", cxxopts::value<std::string>())
    ;
  cli.add_options("Output")
    ("w,outdir", "Path to output directory.",
      cxxopts::value<std::string>()->default_value("./"))
    ("B,dump-binary",
      "Binary Dump mode: write ref. tree in binary format then exit.")
    ("c,bfast",
      "Convert the given fasta file to bfast format needed for running EPA-ng with MPI",
      cxxopts::value<std::string>())
    ("filter-acc-lwr",
      "Accumulated likelihood weight after which further placements are discarded.",
      cxxopts::value<double>())
    ("filter-min-lwr",
      "Minimum likelihood weight below which a placement is discarded.",
      cxxopts::value<double>()->default_value("0.01"))
    ("filter-min",
      "Minimum number of placements per sequence to include in final output.",
      cxxopts::value<unsigned int>()->default_value("1"))
    ("filter-max",
      "Maximum number of placements per sequence to include in final output.",
      cxxopts::value<unsigned int>()->default_value("7"))
    ;
  cli.add_options("Compute")
    // ("O,opt-ref-tree", "Optimize reference tree and model parameters.")
    ("g,dyn-heur",
      "Two-phase heuristic, determination of candidate edges using accumulative threshold. Enabled by default! See --no-heur for disabling it",
      cxxopts::value<double>()->default_value("0.99999")->implicit_value("0.99999"))
    ("G,fix-heur",
      "Two-phase heuristic, determination of candidate edges by specified percentage of total edges.",
      cxxopts::value<double>()->implicit_value("0.1"))
    ("baseball-heur",
      "Baseball heuristic as known from pplacer. strike_box=3,max_strikes=6,max_picthes=40.")
    ("no-heur",
      "Disables heuristic preplacement completely. Overrides all other heuristic flags.")
    ("m,model",
      "Description string of the model to be used. May also be a file containing the parameters, such as a RAxML_info file."
      " --model STRING | FILE "
      "See: https://github.com/amkozlov/raxml-ng/wiki/Input-data#evolutionary-model",
      cxxopts::value<std::string>()->default_value("DNA-GTR-EMPIRICAL"))
    ("raxml-blo",
      "Employ old style of branch length optimization during thorough insertion as opposed to sliding approach. "
      "WARNING: may significantly slow down computation.")
    ("no-repeats",
      "Do NOT employ site repeats optimization. (not recommended, will increase memory footprint without improving runtime or quality) ")
    ("no-pre-mask",
      "Do NOT pre-mask sequences.")
    ;
  cli.add_options("Pipeline")
    ("pipeline",
      "Type of distributed parallelism to use. If specified, pieline mode is used. This mode was built to handle "
      "input reference data that is too large to fit into memory of one node. Less efficient than the standard mode.")
    ("chunk-size",
      "Number of query sequences to be read in at a time. May influence performance.",
      cxxopts::value<unsigned int>()->default_value("5000"))
    #ifdef __OMP
    ("T,threads",
      "Number of threads to use. If 0 is passed as argument, program will run with the maximum number "
      "of threads available.",
      cxxopts::value<unsigned int>()->default_value("0"))
    #endif
    ;

  cli.parse(argc, argv);

  if (cli.count("help") or empty) {
    std::cout << cli.help({"", "Input", "Output", "Compute", "Pipeline"});
    exit_epa();
  }

  if (cli.count("outdir")) {
    work_dir = cli["outdir"].as<std::string>();
  }
  ensure_dir_has_slash(work_dir);

  #ifdef __MPI
  genesis::utils::Logging::log_to_file(work_dir + std::to_string(local_rank) + ".epa_info.log");
  #else
  genesis::utils::Logging::log_to_file(work_dir + "epa_info.log");
  #endif

  LOG_INFO << "Selected: Output dir: " << work_dir;

  if (cli.count("verbose")) {
    LOG_INFO << "Selected: verbose (debug) output";
    genesis::utils::Logging::max_level(genesis::utils::Logging::kDebug2);
  }

  if (cli.count("bfast")) {
    LOG_INFO << "Converting given FASTA file to BFAST format.";
    auto fasta = cli["bfast"].as<std::string>();
    LOG_INFO << "Started " << genesis::utils::current_time();
    auto resultfile = Binary_Fasta::fasta_to_bfast(fasta, work_dir);
    LOG_INFO << "Finished " << genesis::utils::current_time();
    LOG_INFO << "Resulting bfast file was written to: " << resultfile;
    exit_epa();
  }

  // check for valid input combinations
  if (not(
        ( cli.count("tree") and cli.count("ref-msa") )
    or  ( cli.count("binary") and (cli.count("query") or cli.count("ref-msa")) )
    )) {
    LOG_INFO << "Must supply reference tree/msa either directly or as precomputed binary.";
    exit_epa(EXIT_FAILURE);
  }

  if (cli.count("query")) {
    query_file = cli["query"].as<std::string>();
    LOG_INFO << "Selected: Query file: " << query_file;
  }

  if (cli.count("tree")) {
    tree_file = cli["tree"].as<std::string>();
    LOG_INFO << "Selected: Tree file: " << tree_file;
  }

  if (cli.count("ref-msa")) {
    reference_file = cli["ref-msa"].as<std::string>();
    LOG_INFO << "Selected: Reference MSA: " << reference_file;
  }

  if (cli.count("binary")) {
    binary_file = cli["binary"].as<std::string>();
    options.load_binary_mode = true;
    LOG_INFO << "Selected: Binary CLV store: " << binary_file;
  }

  if (cli.count("filter-acc-lwr"))
  {
    options.support_threshold = cli["filter-acc-lwr"].as<double>();
    options.acc_threshold = true;
    LOG_INFO << "Selected: Filtering by accumulated threshold: " << options.support_threshold;
  }

  if (cli.count("filter-min-lwr")) {
    options.support_threshold = cli["filter-min-lwr"].as<double>();
    options.acc_threshold = false;
    LOG_INFO << "Selected: Filtering by minimum threshold: " << options.support_threshold;
  }

  if (cli.count("filter-min")) {
    options.filter_min = cli["filter-min"].as<unsigned int>();
    LOG_INFO << "Selected: Minimum number of placements per query: " << options.filter_min;
  }

  if (cli.count("filter-max")) {
    options.filter_max = cli["filter-max"].as<unsigned int>();
    LOG_INFO << "Selected: Maximum number of placements per query: " << options.filter_max;
  }

  if (options.filter_min > options.filter_max) {
    throw std::runtime_error{"filter-min must not exceed filter-max!"};
  }

  // ensure only one heuristic was selected
  if ( (cli.count("fix-heur")
      + cli.count("dyn-heur")
      + cli.count("baseball-heur")
      + cli.count("no-heur")) <= 1 ) {
    throw std::runtime_error{"Heuristic flags are mutually exclusive! Please select only one."};
  }

  if (cli.count("fix-heur")) {
    options.prescoring_threshold = cli["fix-heur"].as<double>();
    options.prescoring = options.prescoring_by_percentage = true;
    LOG_INFO << "Selected: Prescoring by percentage of branches: " << options.prescoring_threshold;
  }

  if (cli.count("dyn-heur")) {
    options.prescoring_threshold = cli["dyn-heur"].as<double>();
    options.prescoring = true;
    options.prescoring_by_percentage = false;
    LOG_INFO << "Selected: Prescoring by accumulated LWR threshold: " << options.prescoring_threshold;
  }

  if (cli.count("baseball-heur")) {
    options.baseball = true;
    LOG_INFO << "Selected: Prescoring using the baseball heuristic";
  }

  if (cli.count("raxml-blo")) {
    options.sliding_blo = false;
    LOG_INFO << "Selected: On query insertion, optimize branch lengths the way RAxML-EPA did it";
  }

  if (cli.count("no-repeats")) {
    options.repeats = false;
    LOG_INFO << "Selected: Using the non-repeats version of libpll/modules";
  }

  if (cli.count("no-pre-mask")) {
    options.premasking = false;
    LOG_INFO << "Selected: Disabling pre-masking";
  }

  if (cli.count("pipeline")) {
    pipeline = true;
    LOG_INFO << "Selected: Using the pipeline distributed parallel scheme.";
  }

  if (cli.count("no-heur")) {
    options.prescoring = false;
    LOG_INFO << "Selected: Disabling the prescoring heuristics.";
  }

  if (cli.count("dump-binary")) {
    options.dump_binary_mode =  true;
    LOG_INFO << "Selected: Build reference tree and write it out as a binary CLV store (for MPI)";
    LOG_INFO << "\tWARNING: this mode means that no placement will take place in this run";
  }

  if (cli.count("model")) {
    model_desc = cli["model"].as<std::string>();
    if (is_file(model_desc)) {
      LOG_INFO << "Selected: Specified model file: " << model_desc;
      model_desc = parse_model(model_desc);
    } else {
      LOG_INFO << "Selected: Specified model: " << model_desc;
    }
  }

  model = raxml::Model(model_desc);

  LOG_INFO << model;

  if (cli.count("chunk-size")) {
    options.chunk_size = cli["chunk-size"].as<unsigned int>();
    LOG_INFO << "Selected: Reading queries in chunks of: " << options.chunk_size;
  }

  if (cli.count("threads")) {
    options.num_threads = cli["threads"].as<unsigned int>();
    LOG_INFO << "Selected: Using threads: " << options.num_threads;
  }

  } catch (const cxxopts::OptionException& e) {
    std::cout << "error parsing options: " << e.what() << std::endl;
    exit_epa(EXIT_FAILURE);
  }


  //================================================================
  //============    EPA    =========================================
  //================================================================

  banner += "    ______ ____   ___           _   __ ______\n";
  banner += "   / ____// __ \\ /   |         / | / // ____/\n";
  banner += "  / __/  / /_/ // /| | ______ /  |/ // / __  \n";
  banner += " / /___ / ____// ___ |/_____// /|  // /_/ /  \n";
  banner += "/_____//_/    /_/  |_|      /_/ |_/ \\____/   \n";
  banner += " \n";

  LOG_INFO << banner << std::endl;

  LOG_DBG << "Peeking into MSA files and generating masks";

  MSA_Info ref_info;
  if (reference_file.size()) {
    ref_info = make_msa_info(reference_file);
    LOG_DBG << "Reference File:\n" << ref_info;
  }

  MSA_Info qry_info;
  if (query_file.size()) {
    qry_info = make_msa_info(query_file);
    LOG_DBG << "Query File:\n" << qry_info;
  }

  MSA_Info::or_mask(ref_info, qry_info);

  // msa_info.reset_gaps <- --no-pre-mask

  MSA ref_msa;
  if (reference_file.size()) {
    ref_msa = build_MSA_from_file(reference_file, ref_info, options.premasking);
    LOG_DBG << "Reference File size: " << ref_msa.size();
    LOG_DBG << "Reference File width: " << ref_msa.num_sites();
    if (ref_msa.size() == 0 or ref_msa.num_sites() == 0 ) {
      throw std::runtime_error{"Something went wrong reading the reference file."};
    }
  }

  // build the Tree
  Tree tree;
  if (options.load_binary_mode) {
    LOG_INFO << "Loading from binary";
    tree = Tree(binary_file, model, options);
  } else {
    // build the full tree with all possible clv's
    tree = Tree(tree_file, ref_msa, model, options);
  }

  if (not options.dump_binary_mode) {
    if (query_file.size() == 0) {
      throw std::runtime_error{"Must supply query file! Combined MSA files not currently supported, please split them and specify using -s and -q."};
    }
  } else {
    // dump to binary if specified
    LOG_INFO << "Writing to binary";
    std::string dump_file(work_dir + "epa_binary_file");
    dump_to_binary(tree, dump_file);
    exit_epa();
  }

  // start the placement process and write to file
  auto start_place = std::chrono::high_resolution_clock::now();
  if (pipeline) {
    pipeline_place(tree, query_file, qry_info, work_dir, options, invocation);
  } else {
    simple_mpi(tree, query_file, qry_info, work_dir, options, invocation);
  }
  auto end_place = std::chrono::high_resolution_clock::now();
  auto placetime = std::chrono::duration_cast<std::chrono::seconds>(end_place - start_place).count();

  LOG_INFO << "Time spent placing: " << placetime << "s";

  MPI_FINALIZE();

  auto end_all = std::chrono::high_resolution_clock::now();
  auto alltime = std::chrono::duration_cast<std::chrono::seconds>(end_all - start_all).count();

  LOG_INFO << "Elapsed Time: " << alltime << "s";

	return EXIT_SUCCESS;
}
