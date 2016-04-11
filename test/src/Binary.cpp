#include "Epatest.hpp"

#include "src/Tree.hpp"
#include "src/Binary.hpp"
#include "src/Options.hpp"
#include "src/Model.hpp"
#include "src/file_io.hpp"


TEST(Binary, write)
{
  // setup
  auto msa = build_MSA_from_file(env->reference_file);
  Model model;
  Options options;
  Tree tree(env->tree_file, msa, model, options);

  // test
  dump_to_binary(tree, env->binary_file);

}

TEST(Binary, read)
{
  // setup
  auto msa = build_MSA_from_file(env->reference_file);
  Model model;
  double freqs[4] = {0.1,0.2,0.3,0.4};
  model.base_frequencies(freqs, 4);
  Options options;
  Tree original_tree(env->tree_file, msa, model, options);
  dump_to_binary(original_tree, env->binary_file);

  // test
  Tree read_tree(env->binary_file, options);
  auto read_freqs = read_tree.model().base_frequencies();

  for (size_t i = 0; i < 4; i++)
    EXPECT_DOUBLE_EQ(freqs[i], read_freqs[i]);
}
