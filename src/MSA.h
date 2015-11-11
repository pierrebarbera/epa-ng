#include <vector>
#include <string>

class MSA
{
public:
  MSA(const std::string& msa_file);
  ~MSA();
  std::tuple<std::string, std::string> get(int i);

  int num_sites;

private:
  void build_from_file(const std::string& msa_file);

  std::vector<std::string> headers;
  std::vector<std::string> sequences;
};