#pragma once

#include <boost/math/constants/constants.hpp>

#include "fmmtl/meta/integer_sequence.hpp"

/** Computes the ith Chebyshev node of an N-sized quadrature on [-1/2, 1/2] */
template <typename T>
constexpr T chebyshev_node(std::size_t i, std::size_t N) {
  using boost::math::constants::pi;
  return
      (i     == 0  ) ? -0.5 :    // left point
      (i     == N-1) ?  0.5 :    // right point
      (2*i+1 == N  ) ?  0.0 :    // N is odd and i is middle
      (i     <  N/2) ? -std::cos(     i  * pi<T>()/(N-1)) / 2 :  // Force symmetry
      (i     >= N/2) ?  std::cos((N-1-i) * pi<T>()/(N-1)) / 2 :
      T();
}


template <typename T, typename Seq>
struct ChebyshevImpl;

template <typename T, std::size_t... Is>
struct ChebyshevImpl<T,fmmtl::index_sequence<Is...>> {
  static constexpr std::size_t N = sizeof...(Is);
  static constexpr T x[N] = { chebyshev_node<T>(Is, N)... };
};
template <typename T, std::size_t... Is>
constexpr T ChebyshevImpl<T,fmmtl::index_sequence<Is...>>::x[];

/** Precompute N Chebyshev nodes of type T in the range [-1/2, 1/2] */
template <typename T, std::size_t N>
struct Chebyshev : ChebyshevImpl<T,fmmtl::make_index_sequence<N>> {};
