#pragma once
/** @file Vec.hpp
 * @brief A small N-dimensional numerical vector type that works on CPU and GPU.
 */

#include <iostream>
#include <cmath>

#include "fmmtl/config.hpp"

#define for_i for(std::size_t i = 0; i != N; ++i)

/** @class Vec
 * @brief Class representing ND points and vectors.
 */
template <std::size_t N, typename T = double>
struct Vec {
  FMMTL_STATIC_ASSERT(N > 0, "Vec<N,T> needs N >= 1");

  T elem[N];

  typedef T               value_type;
  typedef T&              reference;
  typedef const T&        const_reference;
  typedef T*              iterator;
  typedef const T*        const_iterator;
  typedef std::size_t     size_type;
  typedef std::ptrdiff_t  difference_type;

  // CONSTRUCTORS

  FMMTL_INLINE Vec() {
    for_i elem[i] = value_type();
  }
  // TODO: Force 0-initialization to get POD and trivial semantics?
  //FMMTL_INLINE Vec() = default;
  FMMTL_INLINE explicit Vec(value_type b) {
    for_i elem[i] = b;
  }
  FMMTL_INLINE Vec(value_type b0, value_type b1) {
    FMMTL_STATIC_ASSERT(N >= 2, "Too many arguments to Vec constructor");
    elem[0] = b0; elem[1] = b1;
    for(std::size_t i = 2; i != N; ++i) elem[i] = value_type();
  }
  FMMTL_INLINE Vec(value_type b0, value_type b1, value_type b2) {
    FMMTL_STATIC_ASSERT(N >= 3,  "Too many arguments to Vec constructor");
    elem[0] = b0; elem[1] = b1; elem[2] = b2;
    for(std::size_t i = 3; i != N; ++i) elem[i] = value_type();
  }
  FMMTL_INLINE Vec(value_type b0, value_type b1, value_type b2, value_type b3) {
    FMMTL_STATIC_ASSERT(N >= 4,  "Too many arguments to Vec constructor");
    elem[0] = b0; elem[1] = b1; elem[2] = b2; elem[3] = b3;
    for(std::size_t i = 4; i != N; ++i) elem[i] = value_type();
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
template <std::size_t N, typename T>
inline std::ostream& operator<<(std::ostream& s, const Vec<N,T>& a) {
  for_i s << a[i] << " ";
  return s;
}
/** Read a Vec from an input stream */
template <std::size_t N, typename T>
inline std::istream& operator>>(std::istream& s, Vec<N,T>& a) {
  for_i s >> a[i];
  return s;
}

/** Compute cross product of two 3D Vecs */
template <typename T>
FMMTL_INLINE Vec<3,T> cross(const Vec<3,T>& a, const Vec<3,T>& b) {
  return Vec<3,T>(a[1]*b[2] - a[2]*b[1],
                  a[2]*b[0] - a[0]*b[2],
                  a[0]*b[1] - a[1]*b[0]);
}

// ARITHMETIC OPERATORS

/** Unary negation: Return -@a a */
template <std::size_t N, typename T>
FMMTL_INLINE Vec<N,T> operator-(Vec<N,T> a) {
  for_i a[i] = -a[i];
  return a;
}
/** Unary plus: Return @a a. ("+a" should work if "-a" works.) */
template <std::size_t N, typename T>
FMMTL_INLINE Vec<N,T> operator+(const Vec<N,T>& a) {
  return a;
}
template <std::size_t N, typename T>
FMMTL_INLINE Vec<N,T> operator+(Vec<N,T> a, const Vec<N,T>& b) {
  return a += b;
}
template <std::size_t N, typename T, typename D>
FMMTL_INLINE Vec<N,T> operator+(Vec<N,T> a, const D& b) {
  return a += b;
}
template <std::size_t N, typename T, typename D>
FMMTL_INLINE Vec<N,T> operator+(const D& b, Vec<N,T> a) {
  return a += b;
}
template <std::size_t N, typename T>
FMMTL_INLINE Vec<N,T> operator-(Vec<N,T> a, const Vec<N,T>& b) {
  return a -= b;
}
template <std::size_t N, typename T, typename D>
FMMTL_INLINE Vec<N,T> operator-(Vec<N,T> a, const D& b) {
  return a -= b;
}
template <std::size_t N, typename T, typename D>
FMMTL_INLINE Vec<N,T> operator-(const D& b, const Vec<N,T>& a) {
  return (-a) += b;
}
template <std::size_t N, typename T>
FMMTL_INLINE Vec<N,T> operator*(Vec<N,T> a, const Vec<N,T>& b) {
  return a *= b;
}
template <std::size_t N, typename T, typename D>
FMMTL_INLINE Vec<N,T> operator*(Vec<N,T> a, const D& b) {
  return a *= b;
}
template <std::size_t N, typename T, typename D>
FMMTL_INLINE Vec<N,T> operator*(const D& b, Vec<N,T> a) {
  return a *= b;
}
template <std::size_t N, typename T>
FMMTL_INLINE Vec<N,T> operator/(Vec<N,T> a, const Vec<N,T>& b) {
  return a /= b;
}
template <std::size_t N, typename T, typename D>
FMMTL_INLINE Vec<N,T> operator/(Vec<N,T> a, const D& b) {
  return a /= b;
}

// ELEMENTWISE OPERATORS

template <std::size_t N, typename T>
FMMTL_INLINE Vec<N,T> abs(Vec<N,T> a) {
  using std::abs;
  for_i a[i] = abs(a[i]);
  return a;
}
template <std::size_t N, typename T>
FMMTL_INLINE Vec<N,T> sqrt(Vec<N,T> a) {
  using std::sqrt;
  for_i a[i] = sqrt(a[i]);
  return a;
}

// META OPERATIONS

#include "fmmtl/meta/dimension.hpp"
#include "fmmtl/meta/random.hpp"

namespace fmmtl {

template <std::size_t N, typename T>
struct dimension<Vec<N,T> > {
  const static std::size_t value = N;
};

template <std::size_t N, typename T>
struct random<Vec<N,T> > {
  static Vec<N,T> get(T a, T b) {
    Vec<N,T> v;
    for_i v[i] = fmmtl::random<T>::get(a, b);
    return v;
  }
  static Vec<N,T> get() {
    return get(T(0), T(1));
  }
};

}  // end namespace fmmtl

#undef for_i


#include "fmmtl/numeric/Norm.hpp"
