#pragma once

#include <limits>
#include <type_traits>
#include <iterator>
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
    return std::uniform_int_distribution<T>(a,b)(default_generator);
  }
  static result_type get() {
    return std::uniform_int_distribution<T>()(default_generator);
  }
};

template <typename T>
struct random<T, typename enable_if<is_floating_point<T>::value>::type> {
  typedef T result_type;
  result_type operator()(T a, T b) const { return get(a,b); }
  result_type operator()()         const { return get();    }

  static result_type get(T a, T b) {
    return std::uniform_real_distribution<T>(a,b)(default_generator);
  }
  static result_type get() {
    return std::uniform_real_distribution<T>()(default_generator);
  }
};

template <typename T>
struct random<fmmtl::complex<T> > {
  typedef fmmtl::complex<T> result_type;
  result_type operator()(T a, T b) const { return get(a,b); }
  result_type operator()()         const { return get();    }

  static result_type get(T a, T b) {
    return fmmtl::complex<T>(random<T>::get(a,b), random<T>::get(a,b));
  }
  static result_type get() {
    return fmmtl::complex<T>(random<T>::get(), random<T>::get());
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
    return std::complex<T>(random<T>::get(), random<T>::get());
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
    Vec<N,T> v;
    for (std::size_t i = 0; i != N; ++i)
      v[i] = fmmtl::random<T>::get();
    return v;
  }
};


class random_n {
  template <typename T>
  struct random_iterator : std::iterator_traits<T*> {
    random_iterator(std::size_t i) : i_(i) {}
    T operator*() const { return random<T>::get(); }
    random_iterator& operator++() { ++i_; return *this; }
    int  operator- (const random_iterator& o) const { return i_ - o.i_;  }
    bool operator==(const random_iterator& o) const { return i_ == o.i_; }
    bool operator!=(const random_iterator& o) const { return i_ != o.i_; }
    std::size_t i_;
  };

  std::size_t N_;

 public:
  random_n(const std::size_t& N) : N_(N) {}

  template <typename Container>
  operator Container() const {
    typedef typename Container::value_type value_type;
    return Container(random_iterator<value_type>(0),
                     random_iterator<value_type>(N_));
  }
};


} // end namespace fmmtl
