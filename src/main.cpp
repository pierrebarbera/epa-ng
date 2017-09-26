#include <iostream>
#include <string>
#include <algorithm>
#include <chrono>

#include <cxxopts.hpp>

#include "net/mpihead.hpp"
#include "util/logging.hpp"
#include "util/Options.hpp"
#include "io/Binary_Fasta.hpp"
#include "io/Binary.hpp"
#include "io/file_io.hpp"
#include "tree/Tree.hpp"
#include "core/Model.hpp"
#include "core/place.hpp"
#include "seq/MSA_Stream.hpp"
#include "seq/MSA.hpp"

static void ensure_dir_has_slash(std::string& dir)
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

void exit_epa(int ret=EXIT_SUCCESS)
{
  MPI_FINALIZE();
  std::exit(ret);
}

int main(int argc, char** argv)
{
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
  std::string sequence_type("DNA");
  std::string model_id("GTR");
  std::string sub_matrix("");
  Options options;

  for (int i = 0; i < argc; ++i) {
    invocation += argv[i];
    invocation += " ";
  }

  std::string query_file("");
  std::string work_dir("");
  std::string tree_file("");
  std::string reference_file("");
  std::string binary_file("");

  std::string banner;

  Model model;

  const bool empty = argc == 1;

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
    ("no-repeats",
      "Do NOT employ site repeats optimization. (not recommended, will increase memory footprint without improving runtime or quality) ")
    ("g,dyn-heur",
      "Two-phase heuristic, determination of candidate edges using accumulative threshold. Enabled by default! See --no-heur for disabling it",
      cxxopts::value<double>()->default_value("0.99")->implicit_value("0.99"))
    ("no-heur",
      "Disables heuristic preplacement completely. Overrides all other heuristic flags.")
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

  if (cli.count("verbose")) {
    LOG_INFO << "Selected: verbose (debug) output";
    genesis::utils::Logging::max_level(genesis::utils::Logging::kDebug2);
  }

  if (cli.count("bfast")) {
    LOG_INFO << "Converting given FASTA file to BFAST format.";
    auto fasta = cli["bfast"].as<std::string>();
    LOG_INFO << "Started " << genesis::utils::current_time();
    Binary_Fasta::fasta_to_bfast(fasta);
    LOG_INFO << "Finished " << genesis::utils::current_time();
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

  if (cli.count("outdir")) {
    work_dir = cli["outdir"].as<std::string>();
    LOG_INFO << "Selected: Output dir: " << work_dir;
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

  if (cli.count("fix-heur")) {
    if (cli.count("no-heur")) {
      LOG_WARN << "WARNING: ignoring --fix-heur/-G as it conflicts with: --no-heur";
    } else {
      options.prescoring_threshold = cli["fix-heur"].as<double>();
      options.prescoring = options.prescoring_by_percentage = true;
      LOG_INFO << "Selected: Prescoring by percentage of branches: " << options.prescoring_threshold;
    }
  }

  if (cli.count("dyn-heur")) {
    if (cli.count("no-heur")) {
      LOG_WARN << "WARNING: ignoring --dyn-heur/-g as it conflicts with: --no-heur";
    } else {
      options.prescoring_threshold = cli["dyn-heur"].as<double>();
      if (cli.count("fix-heur"))  {
        LOG_INFO << "Cannot use -G and -g concurrently! Running with -g " << options.prescoring_threshold ;
      }
      options.prescoring = true;
      options.prescoring_by_percentage = false;
      LOG_INFO << "Selected: Prescoring by accumulated LWR threshold: " << options.prescoring_threshold;
    }
  }

  if (cli.count("opt-ref-tree")) {
    options.opt_branches = options.opt_model = true;
    LOG_INFO << "Selected: Optimizing the reference tree branch lengths and model parameters";
    if (options.load_binary_mode) {
      LOG_INFO << "\tWARNING: this option is ignored as a binary CLV store was supplied!";
    }
  }

  if (cli.count("raxml-blo")) {
    options.sliding_blo = false;
    LOG_INFO << "Selected: On query insertion, optimize branch lengths the way RAxML-EPA did it";
  }

  if (cli.count("no-repeats")) {
    options.repeats = false;
    LOG_INFO << "Selected: Using the non-repeats version of libpll/modules";
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
    auto parts = split_by_delimiter(cli["model"].as<std::string>(), "-");
    auto s = parts.size();
    if (s > 3) {
      throw std::runtime_error{"Supplied too many model arguments! Must be 3 or less."};
    } else {
      if (s >=1) {
        sequence_type = parts[0];
      }
      if (s >= 2) {
        model_id = parts[1];
      }
      if (s >= 3) {
        sub_matrix = parts[2];
      }

      LOG_DBG << "Model descriptor: " << sequence_type << " "
      << model_id << " " << sub_matrix << " ";
    }
    LOG_INFO << "Selected: Specified model: " << sequence_type << " "
      << model_id << " " << sub_matrix;
  }

  model = Model(sequence_type, model_id, sub_matrix);

  if (!options.opt_model) {
    if (cli.count("base-freqs")) {
      auto freq_strings = split_by_delimiter(cli["base-freqs"].as<std::string>(), ":");
      std::vector<double> freqs;
      for (auto& s : freq_strings)
        freqs.push_back(std::stod(s));
      model.base_frequencies(freqs);
    }

    if (cli.count("sub-rates")) {
      auto rate_strings = split_by_delimiter(cli["sub-rates"].as<std::string>(), ":");
      std::vector<double> rates;
      for (auto& s : rate_strings)
        rates.push_back(std::stod(s));
      model.substitution_rates(rates);
    }

    if (cli.count("alpha")) {
      model.alpha(cli["alpha"].as<double>());
    }
  }

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

  ensure_dir_has_slash(work_dir);

  #ifdef __MPI
  genesis::utils::Logging::log_to_file(work_dir + std::to_string(local_rank) + ".epa_info.log");
  #else
  genesis::utils::Logging::log_to_file(work_dir + "epa_info.log");
  #endif

  banner += "    ______ ____   ___           _   __ ______\n";
  banner += "   / ____// __ \\ /   |         / | / // ____/\n";
  banner += "  / __/  / /_/ // /| | ______ /  |/ // / __  \n";
  banner += " / /___ / ____// ___ |/_____// /|  // /_/ /  \n";
  banner += "/_____//_/    /_/  |_|      /_/ |_/ \\____/   \n";
  banner += " \n";

  LOG_INFO << banner << std::endl;

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
      throw std::runtime_error{"Combined MSA files not currently supported, please split them and specify using -s and -q."};
      // MSA tmp;
      // split_combined_msa(ref_msa, tmp, tree);
      
    }
  } else {
    // dump to binary if specified
    LOG_INFO << "Writing to binary";
    std::string dump_file(work_dir + "epa_binary_file");
    dump_to_binary(tree, dump_file);
    exit_epa();
  }

  // start the placement process and write to file
  auto start = std::chrono::high_resolution_clock::now();
  // process(tree, queries, work_dir, options, invocation);
  // tmp_pipeline_test(tree, queries, work_dir, options, invocation);
  simple_mpi(tree, query_file, work_dir, options, invocation);
  auto end = std::chrono::high_resolution_clock::now();
  auto runtime = std::chrono::duration_cast<std::chrono::seconds>(end - start).count();

  LOG_INFO << "Time spent placing: " << runtime << "s";

  MPI_FINALIZE();
	return EXIT_SUCCESS;
}
