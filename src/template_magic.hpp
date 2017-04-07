#pragma once

#include <tuple>
#include <utility>

#include "function_traits.hpp"

template < typename T , typename... Ts >
auto head( std::tuple<T,Ts...> t )
{
   return  std::get<0>(t);
}

template < typename T , typename... Ts >
auto tail( std::tuple<T,Ts...> t )
{
   return  std::get<sizeof...(Ts)>(t);
}

template < std::size_t... Ns , typename... Ts >
static auto skip_impl( std::index_sequence<Ns...> , std::tuple<Ts...> t )
{
   return  std::make_tuple( std::get<Ns + 1u>(t)... );
}

template < typename... Ts >
constexpr auto skip( std::tuple<Ts...> t )
{
   return  skip_impl( std::make_index_sequence<sizeof...(Ts) - 1u>() , t );
}

template < std::size_t... Ns , typename... Ts >
static auto rskip_impl( std::index_sequence<Ns...> , std::tuple<Ts...> t )
{
   return  std::make_tuple( std::get<Ns>(t)... );
}

template < typename... Ts >
constexpr auto rskip( std::tuple<Ts...> t )
{
   return  rskip_impl( std::make_index_sequence<sizeof...(Ts) - 1u>() , t );
}

template < class T, class... Ts >
struct trailing_types
{
  using types = typename std::tuple<Ts...>;
};

// base case
template < class I, class... Ts>
struct leading_types_base;

template < std::size_t... I, class... Ts >
struct leading_types_base<std::index_sequence<I...>, Ts...>
{
  using types = typename std::tuple< typename std::tuple_element< I, std::tuple<Ts...> >::type... >;
};

template < class... Ts >
struct leading_types 
  : leading_types_base<std::make_index_sequence<sizeof...(Ts) - 1u>, Ts...>
{
};

template < class...Ts >
struct inner_types
{
  using all_types = typename leading_types< typename trailing_types<Ts...>::types >::types;
};

/**
 * A tuple for_each function, taken from
 * http://stackoverflow.com/a/6894436
 */
template<std::size_t I = 0, class FuncT, class... Tp>
inline typename std::enable_if<I == sizeof...(Tp), void>::type
for_each(std::tuple<Tp...> &, FuncT) // Unused arguments are given no names.
{ }

template<std::size_t I = 0, class FuncT, class... Tp>
inline typename std::enable_if<I < sizeof...(Tp), void>::type
for_each(std::tuple<Tp...>& t, FuncT f)
{
  f(std::get<I>(t));
  for_each<I + 1, FuncT, Tp...>(t, f);
}
