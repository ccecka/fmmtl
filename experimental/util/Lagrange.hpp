#pragma once

#include <array>
#include <algorithm>

#include "Chebyshev.hpp"


/** An interpolater for the tensor range [-1/2, 1/2]^D using Chebyshev nodes
 */
template <std::size_t Q>
struct Lagrange {

  template <typename T>
  static T coeff(const std::size_t& i, const T& p) {
    typedef Chebyshev<T,Q> C;
    const T& xi = C::x[i];

    T result(1);
    auto x_it  = C::x.begin();

    auto x_end = x_it + i;
    for ( ; x_it != x_end; ++x_it)
      result *= (p - *x_it) / (xi - *x_it);

    ++x_it;

    x_end = C::x.end();
    for ( ; x_it != x_end; ++x_it)
      result *= (p - *x_it) / (xi - *x_it);

    return result;
  }

  template <typename MultiIndex, typename Point>
  static typename Point::value_type
  coeff(const MultiIndex& i,
        const Point& p) {
    typedef typename Point::value_type value_type;
    typedef Chebyshev<value_type,Q> C;

    value_type result(1);
    for (std::size_t k = 0; k < i.size(); ++k) {
      for (std::size_t q = 0; q < Q; ++q) {
        if (q != i[k]) {
          result *= (p[k] - C::x[q]) / (C::x[i[k]] - C::x[q]);
        }
      }
    }
    return result;
  }

  template <typename MultiIndex, typename Point>
  static std::vector<typename Point::value_type>
  coeff(const MultiIndex& i,
        const std::vector<Point>& p) {
    typedef typename Point::value_type value_type;
    typedef Chebyshev<value_type,Q> C;

    std::vector<value_type> result(p.size(), 1);

    for (std::size_t k = 0; k != i.size(); ++k) {
      for (std::size_t q = 0; q != Q; ++q) {
        if (q != i[k]) {
          //const value_type xkxq = C::x[i[k]]-C::x[q]; // Don't help compiler!
          for (std::size_t r = 0; r != p.size(); ++r) {
            result[r] *= (p[r][k] - C::x[q]) / (C::x[i[k]] - C::x[q]);
          }
        }
      }
    }
    return result;
  }

  template <typename MultiIndex, std::size_t N, typename Point>
  static std::array<typename Point::value_type,N>
  coeff(const MultiIndex& i,
        const std::array<Point,N>& p) {
    typedef typename Point::value_type value_type;
    typedef Chebyshev<value_type,Q> C;

    std::array<value_type,N> result;
    result.fill(1);

    for (std::size_t k = 0; k < i.size(); ++k) {
      for (std::size_t q = 0; q < Q; ++q) {
        if (q != i[k]) {
          for (std::size_t r = 0; r != N; ++r) {
            result[r] *= (p[r][k] - C::x[q]) / (C::x[i[k]] - C::x[q]);
          }
        }
      }
    }
    return result;
  }

};


/** Represents a Lagrange Matrix L_i(p_j) for multiindex i and point j
 *
 * L{x}_i(p_j) =  if p_j == x_k for some k,  i == j
 *                else,                      prod_{k!=i} (p_j - x_k)/(x_i - x_k)
 *
 * For a multi-index i of dimension d and point p_j od dimension d:
 * L{x}_i(p_j) = L_{x}_i1(p_j1) * ... * L{x}_id(p_jd)
 */
template <std::size_t D, std::size_t Q, typename T = double>
struct LagrangeMatrix {
  using value_type = T;
  using Cheb = Chebyshev<T,Q>;

  //! Precomputed interpolation weights
  static const std::array<T,Q> w;
  //! The points this matrix is defined over
  std::vector<std::array<T,D>> p;
  //! Precomputed   q[i] = prod_{k} (p[i] - x_k)
  //std::vector<std::array<T,D>> q;
  std::vector<std::array<T,D>> q;

  /** Define w[i] = prod_{k != i}(1 / (x_i - x_k))
   * where z_i are the Chebyshev nodes
   */
  static std::array<T,Q> make() {
    std::array<T,Q> w;
    w.fill(T(1));
    for (std::size_t i = 0; i != Q; ++i) {
      for (std::size_t k = 0; k != Q; ++k)
        if (i != k)
          w[i] *= (Cheb::x[i] - Cheb::x[k]);
      w[i] = T(1) / w[i];
    }
    return w;
  }

  template <typename PointIter>
  LagrangeMatrix(PointIter first, PointIter last) {
    p.resize(std::distance(first,last));
    q.resize(p.size());

    auto pit = p.begin();
    auto qit = q.begin();
    for ( ; first != last; ++first, ++pit, ++qit) {
      auto&& point = *first;
      auto&& pi = *pit;

      std::copy(point.begin(), point.end(), pi.begin());

      std::array<T,D>& qi = *qit;
      qi.fill(value_type(1));
      for (std::size_t k = 0; k != Q; ++k) {
        for (std::size_t d = 0; d != D; ++d) {
          qi[d] *= (pi[d] - Cheb::x[k]);
        }
      }
    }
  }

  value_type operator()(const std::size_t i, const std::size_t d,
                        const std::size_t j) const {
    if (std::abs(p[j][d] - Cheb::x[i]) < 1e-200)
      return value_type(1);
    return w[i] * q[j][d] / (p[j][d] - Cheb::x[i]);
  }

  template <typename MultiIndex>
  value_type operator()(const MultiIndex& i, const std::size_t j) const {
    value_type result = operator()(i[0], 0, j);
    for (std::size_t d = 1; d != D; ++d) {
      result *= operator()(i[d], d, j);
    }
    return result;
  }

  // TODO: Define matvec and transposed matvec
};
template <std::size_t D, std::size_t Q, typename T>
const std::array<T,Q> LagrangeMatrix<D,Q,T>::w = LagrangeMatrix<D,Q,T>::make();
