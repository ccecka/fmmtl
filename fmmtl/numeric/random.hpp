#pragma once

#include <limits>

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int_distribution.hpp>
#include <boost/random/uniform_real_distribution.hpp>

#include "fmmtl/numeric/Vec.hpp"
#include "fmmtl/numeric/Complex.hpp"

namespace fmmtl {

static boost::random::mt19937 default_genenerator;

template <typename T>
struct random;

template <>
struct random<double> {
  static double get(double a, double b) {
    boost::random::uniform_real_distribution<double> dist(a, b);
    return dist(default_genenerator);
  }
  static double get() {
    return get(0,1);
  }
};

template <>
struct random<float> {
  static float get(float a, float b) {
    boost::random::uniform_real_distribution<float> dist(a, b);
    return dist(default_genenerator);
  }
  static float get() {
    return get(0,1);
  }
};

template <>
struct random<unsigned> {
  static unsigned get(unsigned a, unsigned b) {
    boost::random::uniform_int_distribution<unsigned> dist(a, b);
    return dist(default_genenerator);
  }
  static unsigned get() {
    return get(0, std::numeric_limits<unsigned>::max());
  }
};

template <>
struct random<int> {
  static int get(int a, int b) {
    boost::random::uniform_int_distribution<int> dist(a, b);
    return dist(default_genenerator);
  }
  static int get() {
    return get(0, std::numeric_limits<int>::max());
  }
};

template <typename T>
struct random<complex<T> > {
  static complex<T> get(T a, T b) {
    return complex<T>(random<T>::get(a,b), random<T>::get(a,b));
  }
  static complex<T> get() {
    return get(T(0), T(1));
  }
};

template <std::size_t N, typename T>
struct random<Vec<N,T> > {
  static Vec<N,T> get(T a, T b) {
    Vec<N,T> v;
    for (std::size_t i = 0; i != N; ++i)
      v[i] = fmmtl::random<T>::get(a, b);
    return v;
  }
  static Vec<N,T> get() {
    return get(T(0), T(1));
  }
};

} // end namespace fmmtl
