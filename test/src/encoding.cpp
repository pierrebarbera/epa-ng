#include "Epatest.hpp"

#include "src/encoding.hpp"

TEST(encoding, 4bit)
{
  FourBit converter;
  const std::string input("AATGCTTCGTAA---NNNATTCBDAVMKWYR");

  auto packed = converter.to_fourbit(input);

  auto unpacked = converter.from_fourbit(packed, input.size());

  EXPECT_STRCASEEQ(input.c_str(), unpacked.c_str());
  // printf("%s\n", input.c_str());
  // printf("%s\n", unpacked.c_str());
}
