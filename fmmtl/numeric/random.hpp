#pragma once

#include <limits>
#include <type_traits>
#include <random>

#include "fmmtl/numeric/Vec.hpp"
#include "fmmtl/numeric/Complex.hpp"

namespace fmmtl {

static std::mt19937 default_generator;

using std::enable_if;
using std::is_integral;
using std::is_floating_point;

template <typename T, class Enable = void>
struct random;

template <typename T>
struct random<T, typename enable_if<is_integral<T>::value>::type> {
  typedef T result_type;
  result_type operator()(T a, T b) const { return get(a,b); }
  result_type operator()()         const { return get();    }

  static result_type get(T a, T b) {
    std::uniform_int_distribution<T> dist(a, b);
    return dist(default_generator);
  }
  static result_type get() {
    return get(T(0), std::numeric_limits<T>::max());
  }
};

template <typename T>
struct random<T, typename enable_if<is_floating_point<T>::value>::type> {
  typedef T result_type;
  result_type operator()(T a, T b) const { return get(a,b); }
  result_type operator()()         const { return get();    }

  static result_type get(T a, T b) {
    std::uniform_real_distribution<T> dist(a, b);
    return dist(default_generator);
  }
  static result_type get() {
    return get(T(0), T(1));
  }
};

template <typename T>
struct random<complex<T> > {
  typedef complex<T> result_type;
  result_type operator()(T a, T b) const { return get(a,b); }
  result_type operator()()         const { return get();    }

  static result_type get(T a, T b) {
    return complex<T>(random<T>::get(a,b), random<T>::get(a,b));
  }
  static result_type get() {
    return get(T(0), T(1));
  }
};

template <typename T>
struct random<std::complex<T> > {
  typedef std::complex<T> result_type;
  result_type operator()(T a, T b) const { return get(a,b); }
  result_type operator()()         const { return get();    }

  static result_type get(T a, T b) {
    return std::complex<T>(random<T>::get(a,b), random<T>::get(a,b));
  }
  static result_type get() {
    return get(T(0), T(1));
  }
};

template <std::size_t N, typename T>
struct random<Vec<N,T> > {
  typedef Vec<N,T> result_type;
  result_type operator()(T a, T b) const { return get(a,b); }
  result_type operator()()         const { return get();    }

  static result_type get(T a, T b) {
    Vec<N,T> v;
    for (std::size_t i = 0; i != N; ++i)
      v[i] = fmmtl::random<T>::get(a, b);
    return v;
  }
  static result_type get() {
    return get(T(0), T(1));
  }
};


class random_n {
  template <typename T>
  struct random_iterator {
    typedef std::input_iterator_tag iterator_category;
    typedef T value_type;
    typedef int difference_type;
    typedef T* pointer;
    typedef T& reference;

    random_iterator(int n = 0) : i_(n) {}

    T operator*() const { return random<T>::get(); }
    random_iterator& operator++() { ++i_; return *this; }
    int operator-(const random_iterator& o) const { return o.i_ - i_; }
    bool operator==(const random_iterator& o) const { return i_ == o.i_; }
    bool operator!=(const random_iterator& o) const { return i_ != o.i_; }

   private:
    std::size_t i_;
  };

  std::size_t N;

 public:
  random_n(const std::size_t& _N) : N(_N) {}

  template <typename Container>
  operator Container() const {
    typedef typename Container::value_type value_type;
    return Container(random_iterator<value_type>(0),
                     random_iterator<value_type>(N));
  }
};


} // end namespace fmmtl
