#include <gtest/gtest.h>
#include <string>
#include <iostream>

#include "Epatest.hpp"
#include "src/logging.hpp"

Epatest* env;
// Log lgr;

int main(int argc, char** argv)
{
  lgr.to_cout(false);
  env = new Epatest();

  // Set data dir using the program path.
  std::string call = argv[0];
  std::size_t found = call.find_last_of("/\\");
  if (found != std::string::npos) {
      env->data_dir = call.substr(0,found) + "/../data/";
  }

  env->tree_file = std::string(env->data_dir);
  env->tree_file += "ref.tre";
  env->reference_file  = std::string(env->data_dir);
  env->reference_file += "aln.fasta";
  env->query_file  = std::string(env->data_dir);
  env->query_file += "query.fasta";
  env->combined_file  = std::string(env->data_dir);
  env->combined_file += "combined.fasta";
  env->binary_file  = std::string(env->data_dir);
  env->binary_file += "persisted.bin";

  ::testing::InitGoogleTest(&argc, argv);

  ::testing::AddGlobalTestEnvironment(env);
  return RUN_ALL_TESTS();
}
