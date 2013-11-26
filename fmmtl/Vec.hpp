#pragma once
/** @file Vec.hpp
 * @brief A wrapper class to provide simple vector operations for
 * primitive types and classes that only require op[].
 */

#include <iostream>
#include <cmath>

#include "config.hpp"

#define for_i for(std::size_t i = 0; i != N; ++i)

/** @class Vec
 * @brief Class representing ND points and vectors.
 */
template <std::size_t N, typename T = double>
struct Vec {
  T elem[N];

  typedef T               value_type;
  typedef T&              reference;
  typedef const T&        const_reference;
  typedef T*              iterator;
  typedef const T*        const_iterator;
  typedef std::size_t     size_type;
  typedef std::ptrdiff_t  difference_type;

  // CONSTRUCTORS

  FMMTL_INLINE Vec() : elem() {
  }
  // TODO: Force 0-initialization to get POD and trivial semantics?
  //FMMTL_INLINE Vec() = default;
  FMMTL_INLINE explicit Vec(value_type b) {
    for_i elem[i] = b;
  }
  FMMTL_INLINE Vec(value_type b0, value_type b1) {
    FMMTL_STATIC_ASSERT(N >= 2);
    elem[0] = b0; elem[1] = b1;
    for(std::size_t i = 2; i != N; ++i) elem[i] = 0;
  }
  FMMTL_INLINE Vec(value_type b0, value_type b1, value_type b2) {
    FMMTL_STATIC_ASSERT(N >= 3);
    elem[0] = b0; elem[1] = b1; elem[2] = b2;
    for(std::size_t i = 3; i != N; ++i) elem[i] = 0;
  }
  FMMTL_INLINE Vec(value_type b0, value_type b1, value_type b2, value_type b3) {
    FMMTL_STATIC_ASSERT(N >= 4);
    elem[0] = b0; elem[1] = b1; elem[2] = b2; elem[3] = b3;
    for(std::size_t i = 4; i != N; ++i) elem[i] = 0;
  }

  // COMPARATORS

  FMMTL_INLINE bool operator==(const Vec& b) const {
    for_i if (elem[i] != b[i]) return false;
    return true;
  }
  FMMTL_INLINE bool operator!=(const Vec& b) const {
    return !(*this == b);
  }

  // MODIFIERS

  /** Add scalar @a b to this Vec */
  template <typename D>
  FMMTL_INLINE Vec& operator+=(const D& b) {
    for_i elem[i] += b;
    return *this;
  }
  /** Subtract scalar @a b from this Vec */
  template <typename D>
  FMMTL_INLINE Vec& operator-=(const D& b) {
    for_i elem[i] -= b;
    return *this;
  }
  /** Scale this Vec up by scalar @a b */
  template <typename D>
  FMMTL_INLINE Vec& operator*=(const D& b) {
    for_i elem[i] *= b;
    return *this;
  }
  /** Scale this Vec down by scalar @a b */
  template <typename D>
  FMMTL_INLINE Vec& operator/=(const D& b) {
    for_i elem[i] /= b;
    return *this;
  }
  /** Add Vec @a b to this Vec */
  FMMTL_INLINE Vec& operator+=(const Vec& b) {
    for_i elem[i] += b[i];
    return *this;
  }
  /** Subtract Vec @a b from this Vec */
  FMMTL_INLINE Vec& operator-=(const Vec& b) {
    for_i elem[i] -= b[i];
    return *this;
  }
  /** Scale this Vec up by factors in @a b */
  FMMTL_INLINE Vec& operator*=(const Vec& b) {
    for_i elem[i] *= b[i];
    return *this;
  }
  /** Scale this Vec down by factors in @a b */
  FMMTL_INLINE Vec& operator/=(const Vec& b) {
    for_i elem[i] /= b[i];
    return *this;
  }
  /** Compute the dot product of this Vec with another Vec */
  FMMTL_INLINE value_type dot(const Vec& b) const {
    value_type d = value_type();
    for_i d += elem[i]*b[i];
    return d;
  }

  // ACCESSORS

  FMMTL_INLINE reference       operator[](size_type i)       { return elem[i]; }
  FMMTL_INLINE const_reference operator[](size_type i) const { return elem[i]; }

  FMMTL_INLINE T*       data()       { return elem; }
  FMMTL_INLINE const T* data() const { return elem; }

  FMMTL_INLINE reference       front()       { return elem[0]; }
  FMMTL_INLINE const_reference front() const { return elem[0]; }
  FMMTL_INLINE reference        back()       { return elem[N-1]; }
  FMMTL_INLINE const_reference  back() const { return elem[N-1]; }

  FMMTL_INLINE static size_type     size() { return N; }
  FMMTL_INLINE static size_type max_size() { return N; }
  FMMTL_INLINE static bool         empty() { return false; }

  // ITERATORS

  FMMTL_INLINE iterator        begin()       { return elem; }
  FMMTL_INLINE const_iterator  begin() const { return elem; }
  FMMTL_INLINE const_iterator cbegin() const { return elem; }

  FMMTL_INLINE iterator          end()       { return elem+N; }
  FMMTL_INLINE const_iterator    end() const { return elem+N; }
  FMMTL_INLINE const_iterator   cend() const { return elem+N; }
};

// OPERATORS

/** Write a Vec to an output stream */
template <std::size_t N, typename P>
inline std::ostream& operator<<(std::ostream& s, const Vec<N,P>& a) {
  for_i s << a[i] << " ";
  return s;
}
/** Read a Vec from an input stream */
template <std::size_t N, typename P>
inline std::istream& operator>>(std::istream& s, Vec<N,P>& a) {
  for_i s >> a[i];
  return s;
}

/** Compute the dot product of two Vecs */
template <std::size_t N, typename P>
FMMTL_INLINE typename Vec<N,P>::value_type dot(const Vec<N,P>& a,
                                               const Vec<N,P>& b) {
  return a.dot(b);
}
/** Compute the dot product of two Vecs */
template <std::size_t N, typename P>
FMMTL_INLINE typename Vec<N,P>::value_type inner_prod(const Vec<N,P>& a,
                                                      const Vec<N,P>& b) {
  return a.dot(b);
}
/** Compute cross product of two 3D Vecs */
template <typename P>
FMMTL_INLINE Vec<3,P> cross(const Vec<3,P>& a, const Vec<3,P>& b) {
  return Vec<3,P>(a[1]*b[2] - a[2]*b[1],
                  a[2]*b[0] - a[0]*b[2],
                  a[0]*b[1] - a[1]*b[0]);
}
/** Compute the squared L2 norm */
template <std::size_t N, typename P>
FMMTL_INLINE typename Vec<N,P>::value_type normSq(const Vec<N,P>& a) {
  return a.dot(a);
}
/** Compute the L2 norm */
template <std::size_t N, typename P>
FMMTL_INLINE typename Vec<N,P>::value_type norm(const Vec<N,P>& a) {
  using std::sqrt;
  return sqrt(normSq(a));
}
/** Compute the L2 norm */
template <std::size_t N, typename P>
FMMTL_INLINE typename Vec<N,P>::value_type norm_2(const Vec<N,P>& a) {
  return norm(a);
}
/** Compute the L1 norm */
template <std::size_t N, typename P>
FMMTL_INLINE typename Vec<N,P>::value_type norm_1(const Vec<N,P>& a) {
  using std::abs;
  typename Vec<N,P>::value_type r = typename Vec<N,P>::value_type();
  for_i r += abs(a[i]);
  return r;
}
/** Compute the L-infinity norm */
template <std::size_t N, typename P>
FMMTL_INLINE typename Vec<N,P>::value_type norm_inf(const Vec<N,P>& a) {
  using std::abs;
  using std::max;
  typename Vec<N,P>::value_type a_max = typename Vec<N,P>::value_type();
  for_i a_max = max(a_max, abs(a[i]));
  return a_max;
}

// ARITHMETIC OPERATORS

/** Unary negation: Return -@a a */
template <std::size_t N, typename P>
FMMTL_INLINE Vec<N,P> operator-(Vec<N,P> a) {
  for_i a[i] = -a[i];
  return a;
}
/** Unary plus: Return @a a. ("+a" should work if "-a" works.) */
template <std::size_t N, typename P>
FMMTL_INLINE Vec<N,P> operator+(const Vec<N,P>& a) {
  return a;
}
template <std::size_t N, typename P>
FMMTL_INLINE Vec<N,P> operator+(Vec<N,P> a, const Vec<N,P>& b) {
  return a += b;
}
template <std::size_t N, typename P, typename D>
FMMTL_INLINE Vec<N,P> operator+(Vec<N,P> a, const D& b) {
  return a += b;
}
template <std::size_t N, typename P, typename D>
FMMTL_INLINE Vec<N,P> operator+(const D& b, Vec<N,P> a) {
  return a += b;
}
template <std::size_t N, typename P>
FMMTL_INLINE Vec<N,P> operator-(Vec<N,P> a, const Vec<N,P>& b) {
  return a -= b;
}
template <std::size_t N, typename P, typename D>
FMMTL_INLINE Vec<N,P> operator-(Vec<N,P> a, const D& b) {
  return a -= b;
}
template <std::size_t N, typename P, typename D>
FMMTL_INLINE Vec<N,P> operator-(const D& b, const Vec<N,P>& a) {
  return (-a) += b;
}
template <std::size_t N, typename P>
FMMTL_INLINE Vec<N,P> operator*(Vec<N,P> a, const Vec<N,P>& b) {
  return a *= b;
}
template <std::size_t N, typename P, typename D>
FMMTL_INLINE Vec<N,P> operator*(Vec<N,P> a, const D& b) {
  return a *= b;
}
template <std::size_t N, typename P, typename D>
FMMTL_INLINE Vec<N,P> operator*(const D& b, Vec<N,P> a) {
  return a *= b;
}
template <std::size_t N, typename P>
FMMTL_INLINE Vec<N,P> operator/(Vec<N,P> a, const Vec<N,P>& b) {
  return a /= b;
}
template <std::size_t N, typename P, typename D>
FMMTL_INLINE Vec<N,P> operator/(Vec<N,P> a, const D& b) {
  return a /= b;
}

// ELEMENTWISE OPERATORS

template <std::size_t N, typename P>
FMMTL_INLINE Vec<N,P> abs(Vec<N,P> a) {
  using std::abs;
  for_i a[i] = abs(a[i]);
  return a;
}
template <std::size_t N, typename P>
FMMTL_INLINE Vec<N,P> sqrt(Vec<N,P> a) {
  using std::sqrt;
  for_i a[i] = sqrt(a[i]);
  return a;
}

// META OPERATIONS

#include "fmmtl/meta/dimension.hpp"
#include "fmmtl/meta/random.hpp"

namespace fmmtl {

template <std::size_t N, typename P>
struct dimension<Vec<N,P> > {
  const static std::size_t value = N;
};

template <std::size_t N, typename P>
struct random<Vec<N,P> > {
  static Vec<N,P> get(P a, P b) {
    Vec<N,P> v;
    for_i v[i] = fmmtl::random<P>::get(a, b);
    return v;
  }
  static Vec<N,P> get() {
    return get(P(0), P(1));
  }
};

}  // end namespace fmmtl


#undef for_i
