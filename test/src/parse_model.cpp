#include "Epatest.hpp"

#include "util/parse_model.hpp"

using namespace std;

TEST(parse_model, parse_model)
{
  cout << parse_model(env->info_file) << endl;
}
