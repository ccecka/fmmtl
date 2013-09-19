/** @file Math.hpp
 * Allows treating POD the same as we treat a Vec
 */

#include "fmmtl/config.hpp"

#include <cmath>

using std::abs;
using std::sqrt;
// DON'T USE std::norm

/** Compute the dot product of two doubles */
FMMTL_INLINE double dot(double a, double b) {
  return a*b;
}
/** Compute the dot product of two floats */
FMMTL_INLINE double dot(float a, float b) {
  return a*b;
}
/** Compute the dot product of two Vecs */
FMMTL_INLINE double inner_prod(double a, double b) {
  return a*b;
}
/** Compute the dot product of two Vecs */
FMMTL_INLINE double inner_prod(float a, float b) {
  return a*b;
}
/** Compute the squared L2 norm */
FMMTL_INLINE double normSq(double a) {
  return a*a;
}
/** Compute the squared L2 norm */
FMMTL_INLINE float normSq(float a) {
  return a*a;
}
/** Compute the L2 norm */
FMMTL_INLINE double norm(double a) {
  return std::abs(a);
}
/** Compute the L2 norm */
FMMTL_INLINE float norm(float a) {
  return std::abs(a);
}
/** Compute the L2 norm */
FMMTL_INLINE double norm_2(double a) {
  return norm(a);
}
/** Compute the L2 norm */
FMMTL_INLINE float norm_2(float a) {
  return norm(a);
}
/** Compute the L1 norm */
FMMTL_INLINE double norm_1(double a) {
  return norm(a);
}
/** Compute the L1 norm */
FMMTL_INLINE float norm_1(float a) {
  return norm(a);
}
/** Compute the L-infinity norm */
FMMTL_INLINE double norm_inf(double a) {
  return norm(a);
}
/** Compute the L-infinity norm */
FMMTL_INLINE float norm_inf(float a) {
  return norm(a);
}
