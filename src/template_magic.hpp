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

/**
 * A function to link two two-tuples of compatible types into one three-tuple
 */

template <class in, class inter, class out>
auto tuple_link(std::tuple<in, inter> a, std::tuple<inter, out> b)
{
  return std::tuple<in, inter, out>(
    std::get<0>(a),
    std::get<1>(a),
    std::get<1>(b)
  );
}

/**
 * building a tuple of the input/output types of a set of Stages
 */
template <class... more_Stages>
struct token_out_types
{
  // using types = std::tuple<>;
};

template <class Stage, class... more_Stages>
struct token_out_types<Stage, more_Stages...>
{
  using types = typename std::tuple<
    typename Stage::out_type, 
    token_out_types<more_Stages...>
  >;
};

template <class Stage, class... more_Stages>
struct token_types_app
{
  using types = typename std::tuple<
    typename Stage::in_type, 
    typename token_out_types<Stage, more_Stages...>::types
  >;
};

template < class I, class StageTuple>
struct token_types_base;

template < std::size_t... I, class StageTuple >
struct token_types_base<std::index_sequence<I...>, StageTuple>
{
  using types = typename std::tuple< 
    typename std::tuple_element<0u, StageTuple>::type::in_type, 
    typename std::tuple_element<I, StageTuple>::type::out_type... 
  >;
};

template < class StageTuple >
struct token_types
  : token_types_base<std::make_index_sequence< std::tuple_size<StageTuple>::value >, StageTuple>
{
};

// Store parameter pack
template<typename... T>
struct pack
{
    static const unsigned int size = sizeof...(T);
};
 
// Get i-th element of parameter pack
template<int n, typename F, typename... T>
struct element_at : public element_at<n-1, T...>
{
};
 
template<typename F, typename... T>
struct element_at<0, F, T...>
{
    typedef F type;
};
 
// Get i-th element of pack
template<int n, typename P>
struct element
{
};
 
template<int n, typename... T>
struct element<n, pack<T...>>
{
    typedef typename element_at<n, T...>::type type;
};

// Concat at left
template<typename a, typename b>
struct tuple_concat_left
{
};
 
template<typename a, typename... b>
struct tuple_concat_left<a, pack<b...>>
{
    typedef pack<a, b...> type;
};