#include "epa.hpp"

#include <stdexcept>
#include <fstream>
#include <chrono>

#include "file_io.hpp"
#include "jplace_util.hpp"
#include "stringify.hpp"
#include "set_manipulators.hpp"
#include "logging.hpp"

using namespace std;

Log lgr;

static void ensure_outdir_has_slash(string& outdir)
{
  if (outdir.length() > 0 && outdir.back() != '/')
    outdir += "/";
}

void epa(const string& tree_file, const string& reference_msa_file, const string& query_msa_file, string& outdir,
                Model model, const Options& options, const string& invocation)
{
	// sanitize input
  file_check(tree_file);
  file_check(reference_msa_file);
  if(query_msa_file.length() > 0)
    file_check(query_msa_file);

  ensure_outdir_has_slash(outdir);

  lgr = Log(outdir + "epa_info.log");

  lgr << "EPA - Evolutionary Placement Algorithm" << endl;
  lgr << "\nInvocation: \n" << invocation << endl;

	// Build the reference tree
  auto ref_msa = build_MSA_from_file(reference_msa_file);

  MSA query_msa;
  if (query_msa_file.size() != 0)
    query_msa = build_MSA_from_file(query_msa_file);

  auto tree = Tree(tree_file, ref_msa, model, options, query_msa);

  // place query sequences
  auto start = chrono::high_resolution_clock::now();
  auto sample = tree.place();
  auto end = chrono::high_resolution_clock::now();
  auto runtime = chrono::duration_cast<chrono::seconds>(end - start).count();

  lgr << "\nTime spent placing: " << runtime << "s" << endl;

  ofstream outfile(outdir + "epa_result.jplace");
  lgr << "\nWriting output to: " << outdir + "epa_result.jplace" << endl;
  outfile << sample_to_jplace_string(sample, invocation, tree.query_msa()) << endl;
  outfile.close();
}

// void epa(const string& binary_file, const string& tree_file, const string &query_msa_file, string& outdir, Options& options)
// {
//   file_check(binary_file);
//   file_check(query_msa_file);
//
//   ensure_outdir_has_slash(outdir);
//
//   auto tree = Tree(bin_file, tree_file, options);
//
// }
