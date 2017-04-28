#include <iostream>
#include <string>
#include <algorithm>
#include <chrono>

#include <cxxopts.hpp>

#include "mpihead.hpp"
#include "logging.hpp"
#include "epa.hpp"

using namespace std;

// Log lgr;

static void ensure_dir_has_slash(string& dir)
{
  if (dir.length() > 0 && dir.back() != '/') {
    dir += "/";
  }
}

static std::vector<std::string> split_by_delimiter(const std::string & text, const std::string delim)
{
  std::vector<std::string> parts;
  size_t start = 0;
  size_t end = 0;

  do
  {
    end = text.find(delim, start);
    end = (end != std::string::npos) ? end : text.length();
    parts.emplace_back(text.substr(start, end - start));
    start = end + delim.length();
  } while (end != std::string::npos and start <= text.length());

  return parts;
}

int main(int argc, char** argv)
{
  genesis::utils::Logging::log_to_stdout();
  genesis::utils::Logging::max_level(genesis::utils::Logging::kDebug1);

  MPI_INIT(&argc, &argv);

  string invocation("");
  string sequence_type("DNA");
  string model_id("GTR");
  string sub_matrix("");
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

  std::string banner;

  Model model;

  try
  {
  cxxopts::Options cli(argv[0], "Massively-Parallel Evolutionary Placement Algorithm");

  cli.add_options()
    ("help", "Display help.")
    ("version", "Display version.")
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
    ("filter-acc-lwr",
      "Accumulated likelihood weight after which further placements are discarded.",
      cxxopts::value<double>()->default_value("0.9999"))
    ("filter-min-lwr",
      "Minimum likelihood weight below which a placement is discarded.",
      cxxopts::value<double>())
    ("filter-min",
      "Minimum number of placements per sequence to include in final output.",
      cxxopts::value<unsigned int>()->default_value("1"))
    ("filter-max",
      "Maximum number of placements per sequence to include in final output.",
      cxxopts::value<unsigned int>())
    ;
  cli.add_options("Compute")
    ("O,opt-ref-tree", "Optimize reference tree and model parameters.")
    ("raxml-blo",
      "Employ old style of branch length optimization during thorough insertion as opposed to sliding approach. "
      "WARNING: may significantly slow down computation.")
    ("g,dyn-heur",
      "Two-phase heuristic, determination of candidate edges using accumulative threshold.",
      cxxopts::value<double>()->implicit_value("0.99"))
    ("G,fix-heur",
      "Two-phase heuristic, determination of candidate edges by specified percentage of total edges.",
      cxxopts::value<double>()->implicit_value("0.1"))
    ("m,model",
      "Description string of the model to be used. Format: "
      "<type>-<symmetries>-<rate/frequency model> Examples: -m DNA-GTR-EMPIRICAL, -m AA-GTR-BLOSUM62",
      cxxopts::value<std::string>()->default_value("DNA-GTR-EMPIRICAL"))
    ("base-freqs",
      "Base frequencies to be used. Must match alphabet size. Overwritten by -O. Example: "
      "--base-freqs 0.2:0.3:0.25:0.25",
      cxxopts::value<std::string>())
    ("sub-rates",
      "Substitution rates to be used. Must correspond to alphabet size (e.g. 6 for DNA). Overwritten by -O. "
      "Order: A-C, A-G, A-T, C-G, C-T, G-T "
      "Example: --sub-rates 0.88:2.0:1.31:0.86:3.48:1.0",
      cxxopts::value<std::string>())
    ("alpha",
      "Alpha parameter to be used. Overwritten by -O. "
      "Example: --alpha 0.634016",
      cxxopts::value<double>())
    ;
  cli.add_options("Pipeline")
    ("chunk-size",
      "Number of query sequences to be read in at a time. May influence performance.",
      cxxopts::value<unsigned int>()->default_value("1000"))
    #ifdef __OMP
    ("T,threads",
      "Number of threads to use. If 0 is passed as argument, program will run with the maximum number "
      "of threads available.",
      cxxopts::value<unsigned int>()->default_value("0"))
    #endif
    ;

  cli.parse(argc, argv);

  if (cli.count("help"))
  {
    std::cout << cli.help({"", "Input", "Output", "Compute", "Pipeline"});
    exit(EXIT_SUCCESS);
  }

  // check for valid input combinations
  if (not(
        ( cli.count("tree") and cli.count("ref-msa") )
    or  ( cli.count("binary") and (cli.count("query") or cli.count("ref-msa")) )
    ))
  {
    LOG_INFO << "Must supply reference tree/msa either directly or as precomputed binary.";
    exit(EXIT_FAILURE);
  }

  if (cli.count("query")) query_file = cli["query"].as<std::string>();
  if (cli.count("outdir")) work_dir = cli["outdir"].as<std::string>();
  if (cli.count("tree")) tree_file = cli["tree"].as<std::string>();
  if (cli.count("ref-msa")) reference_file = cli["ref-msa"].as<std::string>();
  if (cli.count("binary"))
  {
    binary_file = cli["binary"].as<std::string>();
    options.load_binary_mode = true;
  }
  if (cli.count("filter-acc-lwr"))
  {
    options.support_threshold = cli["filter-acc-lwr"].as<double>();
    options.acc_threshold = true;
  }
  if (cli.count("filter-min-lwr")) 
  {
    options.support_threshold = cli["filter-min-lwr"].as<double>();
    options.acc_threshold = false;
  }
  if (cli.count("filter-min")) 
  {
    options.filter_min = cli["filter-min"].as<unsigned int>();
  }
  if (cli.count("filter-max")) 
  {
    options.filter_max = cli["filter-max"].as<unsigned int>();
  }
  if (options.filter_min > options.filter_max)
    throw std::runtime_error{"filter-min must not exceed filter-max!"};
  if (cli.count("fix-heur"))
  {
    options.prescoring_threshold = cli["fix-heur"].as<double>();
    options.prescoring = options.prescoring_by_percentage = true;
  }
  if (cli.count("dyn-heur"))
  {
    options.prescoring_threshold = cli["dyn-heur"].as<double>();
    if (options.prescoring)  {
      LOG_INFO << "Cannot use -G and -g concurrently! Running with -g " << options.prescoring_threshold ;
    }
    options.prescoring = true;
    options.prescoring_by_percentage = false;
  }
  if (cli.count("opt-ref-tree")) options.opt_branches = options.opt_model = true;
  if (cli.count("raxml-blo")) options.sliding_blo = false;
  if (cli.count("dump-binary")) options.dump_binary_mode =  true;
  if (cli.count("model"))
  {
    auto parts = split_by_delimiter(cli["model"].as<std::string>(), "-");
    auto s = parts.size();
    if (s > 3)
      throw runtime_error{"Supplied too many model arguments! Must be 3 or less."};
    else
    {
      if (s >=1)
        sequence_type = parts[0];
      if (s >= 2)
        model_id = parts[1];
      if (s >= 3)
        sub_matrix = parts[2];

      LOG_DBG << "Model descriptor: " << sequence_type << " "
      << model_id << " " << sub_matrix << " ";
    }
  }

  model = Model(sequence_type, model_id, sub_matrix);

  if (!options.opt_model)
  {
    if (cli.count("base-freqs"))
    {
      auto freq_strings = split_by_delimiter(cli["base-freqs"].as<std::string>(), ":");
      std::vector<double> freqs;
      for (auto& s : freq_strings)
        freqs.push_back(std::stod(s));
      model.base_frequencies(freqs);
    }
    if (cli.count("sub-rates"))
    {
      auto rate_strings = split_by_delimiter(cli["sub-rates"].as<std::string>(), ":");
      std::vector<double> rates;
      for (auto& s : rate_strings)
        rates.push_back(std::stod(s));
      model.substitution_rates(rates);
    }
    if (cli.count("alpha")) 
    {
      model.alpha(cli["alpha"].as<double>());
    }
  }

  if (cli.count("chunk-size")) options.chunk_size = cli["chunk-size"].as<unsigned int>();
  if (cli.count("threads")) options.num_threads = cli["threads"].as<unsigned int>();

  } catch (const cxxopts::OptionException& e) {
    std::cout << "error parsing options: " << e.what() << std::endl;
    exit(EXIT_FAILURE);
  }

  //================================================================
  //============    EPA    =========================================
  //================================================================

  ensure_dir_has_slash(work_dir);

  #ifdef __MPI
  int local_rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &local_rank);
  if (local_rank != 0) {
    genesis::utils::Logging::log_to_stdout(false);  
  }
  genesis::utils::Logging::log_to_file(work_dir + std::to_string(local_rank) + ".epa_info.log");
  #else
  genesis::utils::Logging::log_to_file(work_dir + "epa_info.log");
  #endif

  banner +=   " _____        ______ _____        ";
  banner += "\n|  __ \\      |  ____|  __ \\ /\\    ";
  banner += "\n| |__) |_____| |__  | |__) /  \\   ";
  banner += "\n|  ___/______|  __| |  ___/ /\\ \\  ";
  banner += "\n| |          | |____| |  / ____ \\";
  banner += "\n|_|          |______|_| /_/    \\_\\ \n";

  LOG_INFO << banner;

  MSA ref_msa;
  if (reference_file.size()) {
    ref_msa = build_MSA_from_file(reference_file);
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

  // build the query stream
  MSA_Stream queries;
  if (not options.dump_binary_mode) {
    if (query_file.size() != 0) {
      queries = MSA_Stream(query_file, options.chunk_size);
    } else { // attempt to split msa if it is intermingled with (supposed) query sequences
      throw runtime_error{"Combined MSA files not currently supported, please split them and specify using -s and -q."};
      // split_combined_msa(ref_msa, queries, tree);
    }
  }
  // dump to binary if specified
  if (options.dump_binary_mode) {
    LOG_INFO << "Writing to binary";
    string dump_file(work_dir + "epa_binary_file");
    dump_to_binary(tree, dump_file);
    MPI_FINALIZE();
  	return EXIT_SUCCESS;
  }

  // start the placement process and write to file
  auto start = chrono::high_resolution_clock::now();
  // process(tree, queries, work_dir, options, invocation);
  tmp_pipeline_test(tree, queries, work_dir, options, invocation);
  auto end = chrono::high_resolution_clock::now();
  auto runtime = chrono::duration_cast<chrono::seconds>(end - start).count();

  LOG_INFO << "Time spent placing: " << runtime << "s";

  MPI_FINALIZE();
	return EXIT_SUCCESS;
}
