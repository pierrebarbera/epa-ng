#pragma once

#include <string>
#include <fstream>
#include <streambuf>
#include <unordered_map>

#include "util/logging.hpp"

using ssmap_t = std::unordered_map<std::string, std::string>;

static const ssmap_t raxml8_to_raxmlng = {
  {"GTR", "GTR"},
  {"GTRGAMMA", "GTR+G"},
  {"GTRGAMMA", "GTR+G"}
};

static std::string parse(const std::string& full, std::string qry, size_t& pos )
{
  std::string result;
  pos = full.find(qry, pos);
  if (pos != std::string::npos) {
    pos += qry.length();
    auto end = full.find('\n', pos);
    if (end == std::string::npos) {
      throw std::runtime_error{"couldn't find terminating newline?!"};
    } else {
      result = full.substr(pos, end - pos);
    }
    pos = end;
  }
  return result;
}

static std::string translate(const std::string& str, const ssmap_t& map)
{
  std::string result;

  auto token = map.find(str);

  if (token != map.end()) {
    result = token->second;
  } else {
    throw std::runtime_error{std::string()+"couldn't find translation for token: "+str};
  }
  return result;
}

static std::string from_raxml_8(const std::string& file)
{
  std::string model_desc;
  std::ifstream t(file);
  std::string full((std::istreambuf_iterator<char>(t)),
                    std::istreambuf_iterator<char>());
  size_t pos = 0;

  // parse  sub matrix identifier
  model_desc.append(
    translate(
      parse(full, "Substitution Matrix: ", pos),
      raxml8_to_raxmlng
    )
  );

  // parse alpha
  auto alpha = parse(full, "alpha: ", pos);

  // parse rates
  std::string rates;

  rates += "{";
  rates += parse(full, "rate A <-> C: ", pos) + "/";
  rates += parse(full, "rate A <-> G: ", pos) + "/";
  rates += parse(full, "rate A <-> T: ", pos) + "/";
  rates += parse(full, "rate C <-> G: ", pos) + "/";
  rates += parse(full, "rate C <-> T: ", pos) + "/";
  rates += parse(full, "rate G <-> T: ", pos) + "}";

  model_desc.append(rates);

  // parse stationary frequencies
  std::string freqs;

  freqs += "+FU{";
  freqs += parse(full, "freq pi(A): ", pos) + "/";
  freqs += parse(full, "freq pi(C): ", pos) + "/";
  freqs += parse(full, "freq pi(G): ", pos) + "/";
  freqs += parse(full, "freq pi(T): ", pos) + "}";

  model_desc.append(freqs);

  // append alpha
  alpha = "+G4{" + alpha + "}";
  model_desc.append(alpha);

  return model_desc;
}

inline std::string parse_model(const std::string& file)
{
  return from_raxml_8(file);
}