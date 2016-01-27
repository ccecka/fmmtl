#pragma once

#include <array>
#include <algorithm>

#include "fmmtl/meta/functional.hpp"

#include "Chebyshev.hpp"

/** Compute the inverse Barycentric weights at compile-time
 * w_inv[i] = prod_{k!=i} (x[i] - x[k])
 */
/*
template <typename T, std::size_t N>
constexpr T lagrange_weight(std::size_t i) {
  return (i&1 ?
          (i == 0 || i == N-1 ? -2 : -1) :
          (i == 0 || i == N-1 ?  2 :  1));
}
*/

template <typename T, std::size_t N>
constexpr T lagrange_weight(std::size_t i, std::size_t k = 0) {
  using Cheb = Chebyshev<T,N>;
  return
      (k == i) ? lagrange_weight<T,N>(i,k+1) :
      (k <  N) ? lagrange_weight<T,N>(i,k+1) * (Cheb::x[i] - Cheb::x[k]) :
      T(1);
}


template <typename T, typename Seq>
struct LagrangeWeightImpl;

template <typename T, std::size_t... Is>
struct LagrangeWeightImpl<T,fmmtl::index_sequence<Is...>> {
  static constexpr std::size_t N = sizeof...(Is);
  static constexpr T w[N] = { (T(1) / lagrange_weight<T,N>(Is))... };
};
template <typename T, std::size_t... Is>
constexpr T LagrangeWeightImpl<T,fmmtl::index_sequence<Is...>>::w[];

/** Precompute N Chebyshev nodes of type T in the range [-1/2, 1/2] */
template <typename T, std::size_t N>
struct LagrangeWeight : LagrangeWeightImpl<T,fmmtl::make_index_sequence<N>> {};







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

  //! The points this matrix is defined over
  std::vector<std::array<T,D>> p;
  //! Precomputed   q[i] = prod_{k} (p[i] - x_k)
  std::vector<std::array<T,D>> q;

  /** Lagrange Interpolation Matrix constructor
   *
   * @param[in] first,last  A range of elements (points).
   *                 A point is considered to be D elements of type T accessible
   *                 via begin()/end().
   * @param[in] proj  An optional projection operator acting on the element type
   *                 *first. Defaults to identity, but in the non-identity
   *                 case, the result type must obey the point concept above.
   */
  template <class PointIter, class Proj = fmmtl::identity>
  LagrangeMatrix(PointIter first, PointIter last, Proj proj = Proj()) {
    p.resize(std::distance(first,last));
    q.resize(p.size());

    auto pit = p.begin();
    auto qit = q.begin();
    for ( ; first != last; ++first, ++pit, ++qit) {
      auto&& point = proj(*first);

      // Copy the point to p[i]
      std::array<T,D>& pi = *pit;
      std::copy(point.begin(), point.end(), pi.begin());

      // Compute the pre-factors q[i]
      std::array<T,D>& qi = *qit;
      qi.fill(T(1));
      for (std::size_t k = 0; k != Q; ++k) {
        for (std::size_t d = 0; d != D; ++d) {
          qi[d] *= (pi[d] - Chebyshev<T,Q>::x[k]);
        }
      }
    }
  }

  value_type operator()(const std::size_t i, const std::size_t d,
                        const std::size_t j) const {
    if (std::abs(p[j][d] - Chebyshev<T,Q>::x[i]) < 1e-200)
      return value_type(1);
    return LagrangeWeight<T,Q>::w[i] * q[j][d] / (p[j][d] - Chebyshev<T,Q>::x[i]);
  }

  template <class MultiIndex>
  value_type operator()(const MultiIndex& i, const std::size_t j) const {
    value_type result = operator()(i[0], 0, j);
    for (std::size_t d = 1; d != D; ++d) {
      result *= operator()(i[d], d, j);
    }
    return result;
  }
};


/** Type representing the transpose of a LagrangeMatrix. Allows code like
 *    prod(trans(LagrangeM), ...);
 */
template <std::size_t D, std::size_t Q, class T = double>
struct LagrangeMatrixTranspose {
  const LagrangeMatrix<D,Q,T>& L_;
};

template <std::size_t D, std::size_t Q, class T>
LagrangeMatrixTranspose<D,Q,T> trans(const LagrangeMatrix<D,Q,T>& L) {
  return {L};
}

template <std::size_t D, std::size_t Q, class T>
const LagrangeMatrix<D,Q,T>& trans(const LagrangeMatrixTranspose<D,Q,T>& L) {
  return L.L_;
}


#include "TensorIndexGridRange.hpp"   // TEMP?


/** Compute the product of this matrix and a range of values
 * @param[in] in   An iterator to a range of length L.p.size()
 * @param[in] out  An iterator to a range of length Q^D
 */
template <std::size_t D, std::size_t Q, typename T,
          typename InIter, typename OutIter>
void prod(const LagrangeMatrix<D,Q,T>& L,
          InIter in, OutIter out) {
  // Trivial implementation for now

  for (auto&& i : TensorIndexGridRange<D,Q>()) {
    // Accumulate into y
    auto& yi = *out;
    ++out;

    InIter xj_it = in;
    for (std::size_t j = 0; j < L.p.size(); ++j, ++xj_it)
      yi += L(i,j) * (*xj_it);
  }
}

/** Compute the product of this matrix and a range of values
 * @param[in] in   An iterator to a range of length Q^D
 * @param[in] out  An iterator to a range of length L.p.size()
 */
template <std::size_t D, std::size_t Q, typename T,
          typename InIter, typename OutIter>
void prod(const LagrangeMatrixTranspose<D,Q,T>& LT,
          InIter in, OutIter out) {
  // Trivial implementation for now
  const LagrangeMatrix<D,Q,T>& L = LT.L_;

  for (auto&& i : TensorIndexGridRange<D,Q>()) {
    auto& xi = *in;
    ++in;

    OutIter yj_it = out;
    for (std::size_t j = 0; j < L.p.size(); ++j, ++yj_it)
      *yj_it += L(i,j) * xi;
  }
}
