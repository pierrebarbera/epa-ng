#include <gtest/gtest.h>
#include <string>
#include <iostream>

#include "Epatest.hpp"

Epatest* env;

int main(int argc, char** argv)
{
  env = new Epatest();

  // Set data dir using the program path.
  std::string call = argv[0];
  std::size_t found = call.find_last_of("/\\");
  if (found != std::string::npos) {
      env->data_dir = call.substr(0,found) + "/../data/";
  }

  ::testing::InitGoogleTest(&argc, argv);

  ::testing::AddGlobalTestEnvironment(env);
  return RUN_ALL_TESTS();
}
