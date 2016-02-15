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

void epa(string& tree_file, string& reference_msa_file, string& query_msa_file, string& outdir,
                Model model, Options options, string invocation)
{
	// sanitize input
  file_check(tree_file);
  file_check(reference_msa_file);
  if(query_msa_file.length() > 0)
    file_check(query_msa_file);
  if (outdir.length() > 0 && outdir.back() != '/')
    outdir += "/";
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
  auto pquerys = tree.place();
  auto end = chrono::high_resolution_clock::now();
  auto runtime = chrono::duration_cast<chrono::seconds>(end - start).count();

  lgr << "\nTime spent placing: " << runtime << "s" << endl;

  ofstream outfile(outdir + "epa_result.jplace");
  lgr << "\nWriting output to: " << outdir + "epa_result.jplace" << endl;
  outfile << pquery_set_to_jplace_string(pquerys, invocation) << endl;
  outfile.close();
}
