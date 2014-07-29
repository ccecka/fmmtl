#pragma once

#include <array>

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

  template <std::size_t D, typename Point>
  static typename Point::value_type
  coeff(const std::array<std::size_t,D>& i,
        const Point& p) {
    typedef typename Point::value_type value_type;
    typedef Chebyshev<value_type,Q> C;

    value_type result(1);
    for (std::size_t d = 0; d != D; ++d)
      result *= coeff(i[d], p[d]);

    /*
    value_type result(1);
    for (std::size_t k = 0; k < D; ++k) {
      for (std::size_t q = 0; q < Q; ++q) {
        if (q != i[k]) {
          result *= (p[k] - C::x[q]) / (C::x[i[k]] - C::x[q]);
        }
      }
    }
    */
    return result;
  }

  template <std::size_t D, typename Point>
  static std::vector<typename Point::value_type>
  coeff(const std::array<std::size_t,D>& i,
        const std::vector<Point>& p) {
    typedef typename Point::value_type value_type;
    typedef Chebyshev<value_type,Q> C;

    std::vector<value_type> result(p.size());

    for (std::size_t k = 0; k < p.size(); ++k) {
      result[k] = 1;
      for (std::size_t d = 0; d != D; ++d) {
        result[k] *= coeff(i[d], p[k][d]);
      }
    }

    /*
    for (std::size_t k = 0; k < D; ++k) {
      for (std::size_t q = 0; q < Q; ++q) {
        if (q != i[k]) {
          for (std::size_t r = 0; r != p.size(); ++r) {
            result[r] *= (p[r][k] - C::x[q]) / (C::x[i[k]] - C::x[q]);
          }
        }
      }
    }
    */
    return result;
  }

  template <std::size_t D, std::size_t N, typename Point>
  static std::array<typename Point::value_type,N>
  coeff(const std::array<std::size_t,D>& i,
        const std::array<Point,N>& p) {
    typedef typename Point::value_type value_type;
    typedef Chebyshev<value_type,Q> C;

    std::array<value_type,N> result;
    for (std::size_t k = 0; k < p.size(); ++k) {
      result[k] = 1;
      for (std::size_t d = 0; d != D; ++d) {
        result[k] *= coeff(i[d], p[k][d]);
      }
    }

    /*
    for (std::size_t k = 0; k < D; ++k) {
      for (std::size_t q = 0; q < Q; ++q) {
        if (q != i[k]) {
          for (std::size_t r = 0; r != N; ++r) {
            result[r] *= (p[r][k] - C::x[q]) / (C::x[i[k]] - C::x[q]);
          }
        }
      }
    }
    */
    return result;
  }

};
