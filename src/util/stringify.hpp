#pragma once

#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "core/raxml/Model.hpp"

constexpr char NEWL = '\n';

template< class T, class lambda >
static std::string stringify_vec_impl( std::vector< T > const& vec,
                                       lambda toString )
{
  std::ostringstream output;

  size_t i = 0;
  for( auto& e : vec ) {
    output << toString( e );
    if( ++i < vec.size() ) {
      output << ", ";
    }
  }

  return output.str();
}

std::string stringify( raxml::Model& model );

template< class T >
std::string stringify( std::vector< T > const& vec )
{
  return stringify_vec_impl( vec, []( auto const& elem ) {
    return elem;
  } );
}

template< class T >
std::string stringify( const std::vector< std::vector< T > >& vec )
{
  return stringify_vec_impl( vec, []( auto const& elem ) {
    return elem.size();
  } );
}

std::vector< std::string > split_by_delimiter( std::string const& text, std::string const& delim );
