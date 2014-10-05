
#include <array>

#include "fmmtl/meta/integer_range.hpp"


// Greatest common divisor
constexpr long gcd(long x, long y) {
  return ((x%y) == 0) ? y : gcd(y,x%y);
}

constexpr long sum(long n) {
  return n;
}
template <typename... More>
constexpr long sum(long n, More... ns) {
  return n + sum(ns...);
}

constexpr long prod(long n) {
  return n;
}
template <typename... More>
constexpr long prod(long n, More... ns) {
  return n * sum(ns...);
}

// n! = n * (n-1) * (n-2) * ... * 1
constexpr long factorial(long n) {
	return (n <= 1) ? 1 : n * factorial(n-1);
};

// (n)_k = n! / (n-k)! = n * (n-1) * (n-2) * ... * (n-k+1)
constexpr long permutation(long n, long k) {
	return (k == 1) ? n : (n-k+1) * permutation(n,k-1);
};

// n choose k = n! / (k! (n-k)!)
constexpr long combination(long n, long k) {
	return permutation(n,k) / factorial(k);
};

// The index of the monomial x[d]^p in the graded polynomial
template <std::size_t DIM>
constexpr long mono_index(long d, long p) {
  return (d == 0) ?
      combination(p+DIM-1, DIM) :
      combination(p+(DIM-d)-1, DIM-d) + mono_index<DIM>(d-1, p);
}

// The index of the monomial...
constexpr std::size_t poly_index_impl(unsigned, unsigned tot, std::size_t) {
  // cnt == 1 and tot == i
  return tot;
}

template <typename... More>
constexpr std::size_t poly_index_impl(unsigned cnt, unsigned tot,
                                      std::size_t i, More... is) {
  return combination(tot+cnt-1, cnt) + poly_index_impl(cnt-1, tot-i, is...);
}

template <typename... More>
constexpr std::size_t poly_index(More... is) {
  return poly_index_impl(sizeof...(is), sum(is...), is...);
}

// helper class
template <typename T, std::size_t N, std::size_t... I>
std::size_t poly_index(const std::array<T,N>& params, index_sequence<I...>) {
  return poly_index(std::get<I>(params)...);
}

// "return poly_index(params...)"
template <typename T, std::size_t N>
std::size_t poly_index(const std::array<T,N>& params) {
  return poly_index(params, make_index_sequence<N>());
}


/** Structure for data regarding graded lexigraphical ordered polynomial
 * up to order ORDER and with dimension DIM.
 *
 * A graded polynomial for x = {a,b,c} (DIM = 3) and ORDER = 3, is given by
 * out = {1,
 *        a, b, c,
 *        a^2, ab, ac, b^2, bc, c^2,
 *        a^3, a^2b, a^2c, ab^2, abc, ac^2, b^3, b^2c, bc^2, c^3};
 */
template <typename T, std::size_t DIM, std::size_t ORDER>
struct GradedPolynomial {
  static constexpr std::size_t N = combination(ORDER+DIM, DIM);

  // The monomials (in graded lexigraphic order) composing this polynomial
  std::array<T,N> mono;

  GradedPolynomial() {}

  template <typename Point>
  GradedPolynomial(const Point& x, const T& init = T(1)) {
    mono[0] = init;
    auto xai = mono.begin() + 1;

    // Force unroll and compute mono_index at compile time!!

    // For each order
    for (std::size_t p = 1; p <= ORDER; ++p) {
      // For each dimension
      for (std::size_t d = 0; d < DIM; ++d) {

        // Get the start and end
        auto xak = mono.begin() + mono_index<DIM>(d,p-1);
        auto end = mono.begin() + mono_index<DIM>(0,p);
        for ( ; xak < end; ++xai, ++xak)
          *xai = *xak * x[d];
      }
    }
  }

  static constexpr std::size_t size() {
    return N;
  }

  GradedPolynomial& operator+=(const GradedPolynomial& other) {
    for (std::size_t k = 0; k < N; ++k)
      mono[k] += other.mono[k];
    return *this;
  }
};


template <typename T, std::size_t DIM, std::size_t ORDER>
T inner_prod(const GradedPolynomial<T,DIM,ORDER>& a,
             const GradedPolynomial<T,DIM,ORDER>& b) {
  T result = 0;
  for (std::size_t k = 0; k < a.mono.size(); ++k)
    result += a.mono[k] * b.mono[k];
  return result;
}
