#include "Epatest.hpp"

#include "util/parse_model.hpp"

using namespace std;

TEST(parse_model, oldrax_dna)
{
  auto parsed = parse_model(env->data_dir + "modelfiles/rax8_dna");
  std::string actual = "GTR{0.787874/1.821672/1.294006/0.698421/3.034135/1.000000}+FU{0.256465/0.222535/0.308594/0.212406}+G4{0.478218}";
  ASSERT_STREQ(parsed.c_str(), actual.c_str());
}

TEST(parse_model, oldrax_dna_invar)
{
  auto parsed = parse_model(env->data_dir + "modelfiles/rax8_invar");
  std::string actual = "GTR{1.217620/2.720208/1.342850/1.115245/3.313319/1.000000}+FU{0.222438/0.209333/0.259930/0.308299}+IU{0.051355}+G4{0.532224}";
  ASSERT_STREQ(parsed.c_str(), actual.c_str());
}

TEST(parse_model, oldrax_amino_acid)
{
  auto parsed = parse_model(env->data_dir + "modelfiles/rax8_prot");
  std::string actual = "PROTGTR{1.003440/0.000100/2.196009/5.059275/4.560130/5.912979/5.272299/0.699779/0.020243/0.457291/2.383670/3.120039/0.000100/5.403822/22.495228/10.091005/0.000100/0.854092/7.486708/2.797232/0.000100/3.336754/14.440424/0.699469/2.876841/6.187008/0.146659/1.306731/40.784484/2.932864/0.000100/0.546610/4.283376/1.721473/7.351827/0.000100/0.242215/19.432567/2.883241/4.941582/0.000100/10.033898/16.079424/0.091909/0.223185/8.114498/0.000100/0.000100/0.019247/13.487358/4.159476/0.965976/0.907212/0.159272/0.000100/0.519077/19.701531/3.251195/1.851635/0.000100/0.000100/0.019268/0.609953/0.022740/0.773217/1.163731/1.169671/0.000100/0.396903/0.096786/0.785888/0.349871/2.773764/1.856643/0.000100/0.605999/0.000100/0.000100/4.159594/0.000100/11.101237/2.167720/3.372110/6.818292/2.550358/12.631121/1.065164/22.029064/0.057514/1.610956/15.065560/4.472040/0.339767/4.488237/2.790060/2.836442/0.000100/0.000100/0.356764/1.686983/0.000100/0.000100/0.000100/2.102531/0.000100/0.000100/1.573639/1.971493/1.619030/0.000100/0.000100/1.136575/1.577318/0.000100/0.000100/0.346037/1.243377/0.000100/0.593558/5.343049/0.040971/0.000100/0.794127/0.325353/0.594141/0.481006/0.875320/4.863509/1.002356/1.150123/5.988620/0.463365/5.112361/26.224851/1.003883/19.731303/1.177304/17.282051/2.032055/0.000100/0.090724/4.148174/0.000100/0.000100/46.061900/0.088576/34.484594/9.718744/1.963249/0.441787/0.376570/3.032150/1.205978/4.791255/4.618257/0.000100/0.902389/0.765214/3.363404/0.000100/0.000100/0.339229/0.085804/0.000100/0.901420/4.051803/0.000100/2.292349/5.894876/0.252951/0.069199/0.141572/4.900860/33.781619/1.583086/4.310457/2.531992/1.573469/0.000100/0.510724/16.587826/2.264213/0.442391/0.311316/1.683277/0.672040/5.904625/10.245562/1.330998/1.000000}+FU{0.065149/0.054231/0.041608/0.058452/0.023965/0.036826/0.069410/0.052618/0.030732/0.067906/0.092164/0.051878/0.022917/0.045111/0.040413/0.069908/0.072135/0.004367/0.029144/0.071068}+G4{0.563473}";
  ASSERT_STREQ(parsed.c_str(), actual.c_str());
}

TEST(parse_model, ng_dna)
{
  auto parsed = parse_model(env->data_dir + "modelfiles/raxng_dna");
  std::string actual = "GTR{5.56435/19.04/4.65971/2.04432/69.6551/1}+FC+G4m{0.193259}";
  ASSERT_STREQ(parsed.c_str(), actual.c_str());
}

TEST(parse_model, iqtree_dna)
{
  auto parsed = parse_model(env->data_dir + "modelfiles/iqtree_dna_invar");
  std::string actual = "GTR{0.9467/3.2100/1.8644/0.8054/5.5442/1.0000}+FU{0.2415/0.2465/0.3237/0.1884}+IU{0.1257}+G4{0.8042}";
  ASSERT_STREQ(parsed.c_str(), actual.c_str());
}