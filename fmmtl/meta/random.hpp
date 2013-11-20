#pragma once

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int_distribution.hpp>
#include <boost/random/uniform_real_distribution.hpp>

#include <limits>

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

} // end namepsace fmmtl
