#pragma once

#include <string>
#include <fstream>
#include <streambuf>
#include <unordered_map>

#include "util/logging.hpp"

#include "genesis/utils/text/string.hpp"
#include "genesis/utils/core/fs.hpp"

static std::string parse( std::string const& full, std::string const& qry, size_t& pos )
{
  std::string result;
  pos = full.find( qry, pos );
  if (pos != std::string::npos) {
    pos += qry.length();
    auto end = full.find( '\n', pos );
    if ( end == std::string::npos ) {
      throw std::runtime_error{ "couldnt find terminating newline?!" };
    } else {
      result = full.substr( pos, end - pos );
    }
    pos = end;
  } else {
    throw std::invalid_argument{
      std::string( "Couldn't parse model file! (can't find '") +
      qry + "'!)"
    };
  }
  return result;
}

static std::string parse_between( std::string const& full,
                                  std::string const& lhs,
                                  std::string const& rhs,
                                  size_t& pos )
{
  std::string result;
  // get the latter part of what we're looking for
  auto tail = parse( full, lhs, pos);

  // see if the right hand side delimiter is in the remaining string
  size_t rhs_begin = tail.find( rhs );
  if ( rhs_begin == std::string::npos ) {
    throw std::invalid_argument{
      std::string("Couldn't parse model file! (can't find '") + rhs + "'!)"
    };
  } else if ( rhs_begin == 0 ) {
    throw std::runtime_error{
      std::string("Nothing inbetween '") + rhs + "' and '" + lhs + "'?"
    };
  } else {
    result = tail.substr(0, rhs_begin);
  }

  return result;
}

static bool rest_has( std::string const& full, std::string const& qry, size_t pos )
{
  return ( full.find(qry, pos) != std::string::npos );
}

static std::string from_raxml_8( const std::string& file )
{
  std::string model_desc;
  std::ifstream t( file );
  std::string full((std::istreambuf_iterator<char>(t)),
                    std::istreambuf_iterator<char>());
  size_t pos = 0;

  // parse data type
  auto type = parse( full, "DataType: ", pos );
  const bool dna = (type == "DNA");

  std::string sub_mat;

  // parse  sub matrix identifier
  sub_mat += parse( full, "Substitution Matrix: ", pos );

  if ( not dna and sub_mat == "GTR" ) {
    sub_mat = "PROTGTR";
  }

  model_desc.append( sub_mat );

  // check for and parse alpha
  std::string alpha;
  if( rest_has(full, "alpha: ", pos) ) {
    alpha = "+G4{" + parse( full, "alpha: ", pos ) + "}";
  }

  // check for and parse prop invar
  std::string p_inv;
  if( rest_has( full, "invar: ", pos ) ) {
    p_inv = "+IU{" + parse( full, "invar: ", pos ) + "}";
  }

  const std::string chars = dna ? std::string("ACGT") : std::string("ARNDCQEGHILKMFPSTWYV");

  // parse rates
  std::string rates;

  rates += "{";

  for( size_t i = 0; i < chars.length() - 1u; ++i ) {
    for( size_t k = i + 1u; k < chars.length(); ++k ) {
      if( k > 1 ) {
        rates += "/";
      }
      rates += parse( full,
        std::string("rate ") + chars[i] + " <-> " + chars[k] + ": ",
        pos
      );
    }
  }
  rates += "}";

  model_desc.append( rates );

  // parse stationary frequencies
  std::string freqs;

  freqs += "+FU{";
  for( size_t i = 0; i < chars.length(); ++i ) {
    if( i > 0 ) {
      freqs += "/";
    }
    freqs += parse( full,
      std::string("freq pi(") + chars[i] +"): ",
      pos
    );
  }
  freqs += "}";

  model_desc.append(freqs);

  if( not p_inv.empty() ) {
    model_desc.append( p_inv );
  }

  if( not alpha.empty() ) {
    model_desc.append( alpha );
  }

  return model_desc;
}

static std::string from_raxml_ng( const std::string& file )
{
  std::ifstream fs(file);
  std::string line;
  std::getline( fs, line );
  auto parts = genesis::utils::split( line, ",");

  if ( parts.size() <= 1 ) {
    throw std::runtime_error{"Model string in provided file seems wrong."};
  }

  return parts.at(0);
}

static std::string from_iqtree( const std::string& file )
{
  std::string model_desc;
  std::ifstream t(file);
  std::string full((std::istreambuf_iterator<char>(t)),
                    std::istreambuf_iterator<char>());
  size_t pos = 0;

  // parse  sub matrix identifier
  auto iq_model_string = genesis::utils::split( parse( full, "Model of substitution: ", pos), "+" );

  auto sub_mat = iq_model_string[0];

  bool dna = (sub_mat == "GTR");

  model_desc.append( sub_mat );

  const std::string chars = dna ? std::string("ACGT") : std::string("ARNDCQEGHILKMFPSTWYV");

  // parse rates
  std::string rates;

  rates += "{";

  for( size_t i = 0; i < chars.length() - 1u; ++i ) {
    for( size_t k = i + 1u; k < chars.length(); ++k ) {
      if ( k > 1 ) {
        rates += "/";
      }
      rates += parse( full,
        std::string() + chars[i] + "-" + chars[k] + ": ",
        pos
      );
    }
  }
  rates += "}";
  model_desc.append( rates );

  // parse stationary frequencies
  std::string freqs;
  freqs += "+FU{";
  for( size_t i = 0; i < chars.length(); ++i ) {
    if ( i > 0 ) {
        freqs += "/";
      }
    freqs += parse( full,
      std::string("pi(") + chars[i] +") = ",
      pos
    );
  }
  freqs += "}";
  model_desc.append( freqs );

  // check if gamma is used and parse the number of rates
  bool gamma = rest_has(full,"Gamma with ", pos);
  std::string gamma_cats;
  if ( gamma ) {
    gamma_cats = parse_between(full, "Gamma with ", " categories", pos);
  }

  // check if prop invar was set
  if( rest_has( full, "Proportion of invariable sites: ", pos ) ) {
    std::string p_inv;
    p_inv += "+IU{" + parse( full, "Proportion of invariable sites: ", pos ) + "}";
    model_desc.append( p_inv );
  }

  // parse alpha
  if ( gamma ) {
    auto alpha = parse( full, "Gamma shape alpha: ", pos );
    alpha = "+G" + gamma_cats + "{" + alpha + "}";
    model_desc.append( alpha );
  }

  return model_desc;
}

static bool contains( const std::string& file, const std::string& needle )
{
  auto full_file = genesis::utils::file_read( file );

  return full_file.find( needle ) != std::string::npos;
}

static bool begins_with( const std::string& file, const std::string& needle )
{
  std::ifstream fs( file );
  std::string line;
  std::getline( fs, line );

  return genesis::utils::starts_with( line, needle );
}

enum class ModelfileType {
  kRaxml8,
  kRaxmlng,
  kIqtree,
  kUnknown
};

static ModelfileType guess( const std::string& file )
{
  if ( begins_with( file, "IQ-TREE " ) ) {
    LOG_DBG << "Detected IQ-TREE infofile";
    return ModelfileType::kIqtree;
  } else if ( contains(file, "This is RAxML version 8.") ) {
    LOG_DBG << "Detected RAxML 8 infofile";
    return ModelfileType::kRaxml8;
  } else {
    return ModelfileType::kRaxmlng;
  }
}

std::string parse_model( const std::string& file )
{
  auto filetype = guess( file );

  std::string model_desc;

  switch ( filetype ) {
    case ModelfileType::kRaxmlng:
      model_desc = from_raxml_ng( file );
      break;
    case ModelfileType::kRaxml8:
      model_desc = from_raxml_8( file );
      break;
    case ModelfileType::kIqtree:
      model_desc = from_iqtree( file );
      break;
    case ModelfileType::kUnknown:
      throw std::runtime_error{
        "Could not detect type of model file! Please supply either a RAxML_info file generated "
        " with RAxML v8 (-f e option) or the .bestModel file generated by raxml-ng!"
      };
      break;
    default:
      throw std::runtime_error{"Filetype guess return makes no sense or missing enum value"};
      break;
  }

  return model_desc;
}