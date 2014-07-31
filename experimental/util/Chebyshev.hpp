#pragma once

#include <array>

#include <boost/math/special_functions/cos_pi.hpp>

#include "fmmtl/numeric/Vec.hpp"

/** Precompute N Chebyshev nodes of type T in the range [-1/2, 1/2] */
template <typename T, std::size_t N>
struct Chebyshev {
  static_assert(N > 1, "More than 1 please");
  //! Nodes of the Chebyshev grid
  static const std::array<T,N> x;

  /** Define the nodes of the Chebyshev grid */
  static std::array<T,N> make() {
    std::array<T,N> x;
    for (std::size_t k = 0; k != N; ++k)
      x[k] = boost::math::cos_pi(T(k)/(N-1)) / T(2);
    return x;
  }
};
template <typename T, std::size_t N>
const std::array<T,N> Chebyshev<T,N>::x = Chebyshev<T,N>::make();



/** A lazy Chebyshev grid with custom center and width */
template <std::size_t Q, std::size_t DIM = 1, typename T = double>
struct ChebyshevGrid {
  typedef T value_type;

  typedef std::array<T,Q> container_type;
  typedef typename container_type::const_iterator const_iterator;
  typedef typename container_type::iterator iterator;

  container_type nodes_;   // XXX: Precompute and/or static+const initialization

  /** Default constructor
   * Define Chebyshev grid with Q points on [-1/2, 1/2]
   */
  ChebyshevGrid() : nodes_(Chebyshev<T,Q>::x) {
  }
  /** Shifted constructor
   * Define Chebyshev grid with Q points on [center-width/2,center+width/2]
   */
  ChebyshevGrid(const T& center, const T& width) {
    for (std::size_t k = 0; k != Q; ++k)
      nodes_[k] = center + Chebyshev<T,Q>::x[k] * width;
  }

  constexpr std::size_t size() const {
    return Q;
  }

  const T& operator[](unsigned i) const {
    return nodes_[i];
  }
  const T& operator()(unsigned i) const {
    return nodes_[i];
  }

  const_iterator begin() const { return nodes_.begin(); }
  const_iterator end()   const { return nodes_.end();   }
};
