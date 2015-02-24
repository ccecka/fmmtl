#pragma once
/** @file math.hpp
 * A collection of constexpr math operations.
 */

namespace fmmtl {

// Greatest common divisor
template <typename T>
constexpr T gcd(const T& x, const T& y) {
  return ((x%y) == 0) ? y : gcd(y,x%y);
}

// Sum
template <typename T>
constexpr const T& sum(const T& n) {
  return n;
}
template <typename T, typename... More>
constexpr T sum(const T& n, More... ns) {
  return n + sum(ns...);
}

// Product
template <typename T>
constexpr const T& product(const T& n) {
  return n;
}
template <typename T, typename... More>
constexpr T product(const T& n, More... ns) {
  return n * product(ns...);
}

// Factorial
// n! = n * (n-1) * (n-2) * ... * 1
template <typename T>
constexpr double factorial(const T& n) {
	return (n <= 1) ? 1.0d : n * factorial(n-1);
};

// Combinatorial Permutation
// (n)_k = n! / (n-k)! = n * (n-1) * (n-2) * ... * (n-k+1)
template <typename N, typename K>
constexpr double permutation(const N& n, const K& k) {
	return (k <= 0) ? 1.0d : n * permutation(n-1,k-1);
};

// Combinatorial Combination -- Binomial Coefficient
// n choose k = n! / (k! (n-k)!)
template <typename N, typename K>
constexpr double combination_impl(const N& n, const K& k) {
  return (k <= 0) ? 1.0d : n * combination_impl(n-1,k-1) / k;
}
template <typename N, typename K>
constexpr double combination(const N& n, const K& k) {
  return (n/2 > k) ? combination_impl(n,n-k) : combination_impl(n,k);
};

} // end namespace fmmtl
