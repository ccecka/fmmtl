#pragma once

#include <numeric>
#include <iostream>

#include "fmmtl/meta/integer_sequence.hpp"
#include "fmmtl/meta/for_each.hpp"
#include "fmmtl/meta/math.hpp"
#include "fmmtl/meta/functional.hpp"

namespace fmmtl {

constexpr std::size_t
monomial_index_impl(std::size_t d, std::size_t i) {
  return combination(i+d-1, d);
}
template <typename... More>
constexpr std::size_t
monomial_index_impl(std::size_t d, std::size_t i, std::size_t j, More... is) {
  return combination(i+d-1, d) + monomial_index_impl(d+1, i+j, is...);
}
template <typename... More>
constexpr std::size_t
monomial_index(More... is) {
  return monomial_index_impl(1, is...);
}

constexpr std::size_t last_nz_(std::size_t c) {
  return c;
}
template <typename... More>
constexpr std::size_t last_nz_(std::size_t c, std::size_t i, More... is) {
  return (i!=0) ? last_nz_(1+sizeof...(is), is...) : last_nz_(c, is...);
}
template <typename... More>
constexpr std::size_t last_non_zero(More... is) {
  return last_nz_(0, is...);
}

/** A MultiIndex represents a sequence of integers (n0,n1,...)
 * TODO: Use 8-bit/16-bit words when possible?
 */
template <std::size_t... Is>
struct MultiIndex : index_sequence<Is...> {
  static constexpr std::size_t I   = monomial_index(Is...);
  static constexpr double      F   = product(factorial(Is)...);
  static constexpr std::size_t D   = sizeof...(Is);
  static constexpr std::size_t R   = sum(Is...);
  static constexpr std::size_t TD  = D -  last_non_zero(Is...);
  static constexpr std::size_t TDV = at<TD>(Is...);
};

template <std::size_t... Is, std::size_t... Js>
MultiIndex<(Is+Js)...>
operator+(MultiIndex<Is...>, MultiIndex<Js...>) {
  return {};
}

template <std::size_t... Is, std::size_t... Js>
MultiIndex<(Is > Js ? Is-Js : 0)...>
operator-(MultiIndex<Is...>, MultiIndex<Js...>) {
  return {};
}

template <std::size_t... Is, std::size_t... Js>
constexpr bool
operator<=(MultiIndex<Is...>, MultiIndex<Js...>) {
  return std::integral_constant<bool, product((Is <= Js)...)>::value;
}

template <std::size_t N, std::size_t... Is>
constexpr std::size_t get(MultiIndex<Is...>) {
  return intseq_element_t<N,index_sequence<Is...>>::value;
}

std::ostream& operator<<(std::ostream& s, fmmtl::index_sequence<>) {
  return s;
}
template <std::size_t i, std::size_t... Is>
std::ostream& operator<<(std::ostream& s, fmmtl::index_sequence<i, Is...>) {
  return s << i << " " << fmmtl::index_sequence<Is...>{};
}

/////////////////
// as_multiindex
template <typename S>
struct as_multiindex;

template <std::size_t... Is>
struct as_multiindex<index_sequence<Is...>> {
  using type = MultiIndex<Is...>;
};

template <typename S>
using as_multiindex_t = typename as_multiindex<S>::type;

template <std::size_t DIM, std::size_t D>
using MultiIndexUnit = as_multiindex_t<intseq_cat_t<idxseq_repeat_t<D>,
                                                    index_sequence<1>,
                                                    idxseq_repeat_t<DIM-D-1>>>;

/** Represents a sequence of multiindexes in graded reverse lexigraphical order.
 * The sequence is parameterized:
 *   ORDER: The maximum order of the multiindex -- sum(I...).
 *   DIM:   The dimension of each multiindex    -- sizeof...(I).
 *
 * Example: For ORDER = 4 and DIM = 3, the reverse lexigraphical order is
 *
 *      #  monomial   expon  step
 *     --  ---------  -----  -----
 *      1     1       0 0 0   ---
 *
 *      2        z    0 0 1  (1)*z
 *      3     y       0 1 0  (1)*y
 *      4  x          1 0 0  (1)*x
 *
 *      5        z^2  0 0 2  (2)*z
 *      6     y  z    0 1 1  (3)*z
 *      7  x     z    1 0 1  (4)*z
 *      8     y^2     0 2 0  (3)*y
 *      9  x  y       1 1 0  (4)*y
 *     10  x^2        2 0 0  (4)*x
 *
 *     11        z^3  0 0 3  (5)*z
 *     12     y  z^2  0 1 2  (6)*z
 *     13  x     z^2  1 0 2  (7)*z
 *     14     y^2z    0 2 1  (8)*z
 *     15  x  y  z    1 1 1  (9)*z
 *     16  x^2   z    2 0 1 (10)*z
 *     17     y^3     0 3 0  (8)*y
 *     18  x  y^2     1 2 0  (9)*y
 *     19  x^2y       2 1 0 (10)*y
 *     20  x^3        3 0 0 (10)*x
 *
 *     21        z^4  0 0 4 (11)*z
 *     22     y  z^3  0 1 3 (12)*z
 *     23  x     z^3  1 0 3 (13)*z
 *     24     y^2z^2  0 2 2 (14)*z
 *     25  x  y  z^2  1 1 2 (15)*z
 *     26  x^2   z^2  2 0 2 (16)*z
 *     27     y^3z^1  0 3 1 (17)*z
 *     28  x  y^2z    1 2 1 (18)*z
 *     29  x^2y  z    2 1 1 (19)*z
 *     30  x^3   z    3 0 1 (20)*z
 *     31     y^4     0 4 0 (17)*y
 *     32  x  y^3     1 3 0 (18)*y
 *     33  x^2y^2     2 2 0 (19)*y
 *     34  x^3y       3 1 0 (20)*y
 *     35  x^4        4 0 0 (20)*x
 *
 * Usage:
 * The GradedMonomialSequence may be used to generate a static for-loop as:
 *   fmmtl::for_each(GradedMonomialSequence<4,3>{}, my_function);
 * which calls my_function with a MultiIndex<...> in the order defined above.
 */
template <std::size_t ORDER, std::size_t DIM>
struct GradedMonomialSequence {
  static_assert(DIM > 0, "GradedPolynomial DIM must be > 0");

  template <std::size_t... Is>
  static constexpr std::size_t index(MultiIndex<Is...>) {
    return monomial_index(Is...);
  }
};

template <typename Seq>
struct GMS_Iterator {};

///////////////////////////////////////
// Static Forward Sequence Utilities //
///////////////////////////////////////

// size GradedMonomialSequence
template <std::size_t ORDER, std::size_t DIM>
struct size<GradedMonomialSequence<ORDER,DIM>>
    : std::integral_constant<std::size_t,
                             std::size_t(combination(ORDER+DIM, DIM))> {
};

// begin GradedMonomialSequence  <0, 0, 0, ...>
template <std::size_t ORDER, std::size_t DIM>
struct begin<GradedMonomialSequence<ORDER,DIM>> {
  using type = GMS_Iterator<idxseq_repeat_t<DIM,0>>;
};

// end GradedMonomialSequence    <0, 0, ..., 0, ORDER+1>
template <std::size_t ORDER, std::size_t DIM>
struct end<GradedMonomialSequence<ORDER,DIM>> {
  using type = GMS_Iterator<intseq_push_back_t<idxseq_repeat_t<DIM-1,0>,
                                               index_t<ORDER+1>>>;
};

// deref GMS_Iterator
template <typename Seq>
struct deref<GMS_Iterator<Seq>> {
  using type = as_multiindex_t<Seq>;
};

// next GMS_Iterator
template <std::size_t... I>
struct next_impl;

template <std::size_t i>
struct next_impl<i> {
  using type = index_sequence<i+1>;
};

template <std::size_t i, std::size_t... I>
struct next_impl<i,0,I...> {
  using type = intseq_cat_t<index_sequence<0>, typename next_impl<i,I...>::type>;
};

template <std::size_t i, std::size_t j, std::size_t... I>
struct next_impl<i,j,I...> {
  using type = index_sequence<i+1,j-1,I...>;
};

template <std::size_t i, std::size_t... I>
struct next<GMS_Iterator<index_sequence<i,I...>>> {
  using type = GMS_Iterator<typename next_impl<i,I...>::type>;
};

/*********************/
/** GradedPolynomial */
/*********************/
// TODO: Need C++14 template lambdas...

using std::get;


template <typename ArrayT, typename ArrayX>
struct _pow {
  template <typename N>      // N is a MultiIndex<I...>
  constexpr void operator()(N) {
    using P = decltype(N{} - MultiIndexUnit<N::D,N::TD>{});
    get<N::I>(t) = get<P::I>(t) * get<N::TD>(x) / N::TDV;
  }
  ArrayT& t;
  const ArrayX& x;
};

template <typename ArrayA, typename ArrayB>
struct _plus_eq {
  template <typename N>
  constexpr void operator()(N) {
    get<N::I>(a) += get<N::I>(b);
  }
  ArrayA& a;
  const ArrayB& b;
};

// Accumulate N-K (inner loops)
template <typename N, typename ArrayT, typename ArrayR, typename ArrayM>
struct _acc_nmk {
  template <typename K>
  constexpr typename std::enable_if<K{} <= N{}>::type
  operator()(K) {
    using NmK = decltype(N{} - K{});
    get<N::I>(t) += get<K::I>(r) * get<NmK::I>(m);
  }
  template <typename K>
  constexpr typename std::enable_if<!(K{} <= N{})>::type
  operator()(K) {
  }
  ArrayT& t;
  const ArrayR& r;
  const ArrayM& m;
};

// Accumulate N+K (inner loops)
template <typename N, typename ArrayT, typename ArrayR, typename ArrayM>
struct _acc_npk {
  template <typename K>
  constexpr void operator()(K) {
    using NpK = decltype(N{} + K{});
    get<N::I>(t) += get<K::I>(r) * get<NpK::I>(m);
  }
  ArrayT& t;
  const ArrayR& r;
  const ArrayM& m;
};



template <typename ArrayA, typename ArrayR, typename ArrayM>
struct _sum_to_n_mag {
  template <typename N>
  constexpr void operator()(N) {
    using Seq = GradedMonomialSequence<N::R, N::D>;
    using F   = _acc_nmk<N,ArrayA,ArrayR,ArrayM>;
    for_each<begin_t<Seq>, end_t<Seq>>(F{a,r,m});
  }
  ArrayA& a;
  const ArrayR& r;
  const ArrayM& m;
};

template <typename ArrayA, typename ArrayR, typename ArrayM>
_sum_to_n_mag<ArrayA, ArrayR, ArrayM>
sum_to_n_mag(ArrayA& a, const ArrayR& r, const ArrayM& m) {
  return {a,r,m};
}

template <std::size_t ORDER,
          typename ArrayA, typename ArrayR = ArrayA, typename ArrayM = ArrayA>
struct _sum_to_npk_mag {
  template <typename N>
  constexpr void operator()(N) {
    using Seq = GradedMonomialSequence<ORDER - N::R, N::D>;
    using F   = _acc_npk<N,ArrayA,ArrayR,ArrayM>;
    for_each<begin_t<Seq>, end_t<Seq>>(F{a,r,m});
  }
  ArrayA& a;
  const ArrayR& r;
  const ArrayM& m;
};

template <std::size_t ORDER, typename ArrayA, typename ArrayR, typename ArrayM>
_sum_to_npk_mag<ORDER, ArrayA, ArrayR, ArrayM>
sum_to_npk_mag(ArrayA& a, const ArrayR& r, const ArrayM& m) {
  return {a,r,m};
}

template <std::size_t P, typename ArrayA, typename ArrayR, typename ArrayM>
struct _sum_to_k_mag {
  template <typename N>
  void operator()(N) {
    using Seq = GradedMonomialSequence<P, N::D>;
    using F   = _acc_npk<N,ArrayA,ArrayR,ArrayM>;
    for_each<begin_t<Seq>, end_t<Seq>>(F{a,r,m});
  }

  ArrayA& a;
  const ArrayR& r;
  const ArrayM& m;
};

template <std::size_t ORDER, typename ArrayA, typename ArrayR, typename ArrayM>
_sum_to_k_mag<ORDER, ArrayA, ArrayR, ArrayM>
sum_to_k_mag(ArrayA& a, const ArrayR& r, const ArrayM& m) {
  return {a,r,m};
}


/** Quickie expression template for nice syntax */
template <std::size_t ORDER>
struct multiindex {};

struct gp_mult {};
struct gp_div  {};

template <typename A1, typename Op, typename A2>
struct gp_et {
  const A1& a1;
  const A2& a2;
};

template <typename ArrayLike, std::size_t ORDER>
struct pow_et {
  const ArrayLike& r;

  template <typename T>
  gp_et<pow_et, gp_mult, T> operator*(const T& t) const { return {*this, t}; }
  template <typename T>
  gp_et<pow_et, gp_div, T>  operator/(const T& t) const { return {*this, t}; }
};

template <typename ArrayLike, std::size_t ORDER>
pow_et<ArrayLike,ORDER> pow(const ArrayLike& r, multiindex<ORDER>) {
  return {r};
}



template <typename T, std::size_t D, std::size_t P>
struct GradedPolynomial {
  static constexpr std::size_t ORDER = P;
  static constexpr std::size_t DIM   = D;
  using sequence = GradedMonomialSequence<ORDER,DIM>;
  static constexpr std::size_t N = size<sequence>::value;

  // The monomials (in reverse graded lexigraphic order) of this polynomial
  using MonoArray = std::array<T,N>;
  MonoArray m_;

  GradedPolynomial() = default;

  // Generate the power polynomial
  template <typename ArrayLike>
  GradedPolynomial(const ArrayLike& a, const T& a0) {
    // Set the <0,0,...,0> monomial to the constant
    get<0>(m_) = a0;
    // Skip the <0,0,...,0> monomial...
    using second = next_t<begin_t<sequence> >;
    using last   = end_t<sequence>;
    for_each(second{}, last{}, _pow<MonoArray,ArrayLike>{m_, a});
  }

  // Fill with a constant value
  void fill(const T& value) {
    m_.fill(value);
  }

  static constexpr std::size_t size() noexcept { return N; }

  typename MonoArray::iterator       begin()       { return m_.begin(); }
  typename MonoArray::const_iterator begin() const { return m_.begin(); }

  typename MonoArray::iterator       end()       { return m_.end(); }
  typename MonoArray::const_iterator end() const { return m_.end(); }

        T& operator[](std::size_t i)       { return m_[i]; }
  const T& operator[](std::size_t i) const { return m_[i]; }

  /** Compute this[n] += 1/n! x^n * c for all multiindices n with |n| <= ORDER.
   */
  template <typename Array>
  GradedPolynomial&
  operator+=(const gp_et<pow_et<Array,ORDER>,gp_mult,T>& x) {
    static_assert(fmmtl::dimension<Array>::value == DIM, "Array::size() != DIM");

    // Precompute the multiindex power: c x^m / m!
    GradedPolynomial t(x.a1.r, x.a2);

    // Accumulate into m_, TODO: Fuse with recurrence?
    using first = begin_t<sequence>;
    using last  = end_t<sequence>;
    for_each<first, last>(_plus_eq<MonoArray,MonoArray>{m_, t.m_});
    return *this;
  }

  /** Compute this[n] += sum_k 1/k! x^{k} other[n-k]
   */
  template <typename Array>
  GradedPolynomial&
  operator+=(const gp_et<pow_et<Array,ORDER>,gp_mult,GradedPolynomial>& x) {
    static_assert(fmmtl::dimension<Array>::value == DIM, "Array::size() != DIM");

    // Precompute the multiindex power: x^m / m!
    GradedPolynomial t(x.a1.r, T{1});

    // Now compute the matvec... TODO: Fuse with recurrence?
    using first = begin_t<sequence>;
    using last  = end_t<sequence>;
    for_each<first, last>(sum_to_n_mag(m_, t.m_, x.a2.m_));
    return *this;
  }

  /** Compute this[n] += sum_k 1/k! x^{k} other[n+k]
   */
  template <typename Array>
  GradedPolynomial&
  operator+=(const gp_et<pow_et<Array,ORDER>,gp_div,GradedPolynomial>& x) {
    static_assert(fmmtl::dimension<Array>::value == DIM, "Array::size() != DIM");

    // Precompute the multiindex power: x^m / m!
    GradedPolynomial t(x.a1.r, T{1});

    // Now compute the matvec... TODO: Fuse with recurrence?
    using first = begin_t<sequence>;
    using last  = end_t<sequence>;
    for_each<first, last>(sum_to_npk_mag<ORDER>(m_, t.m_, x.a2.m_));
    return *this;
  }

#if 0
  /** Compute this[n] += sum_k gp1[k] * gp2[n-k]
   */
  template <typename T1, std::size_t P1,
            typename T2, std::size_t P2>
  GradedPolynomial&
  operator+=(const gp_et<GradedPolynomial<T1,D,P1>, gp_mult, GradedPolynomial<T2,D,P2> >& x) {

  }

  /** Compute this[n] += sum_k gp1[k] * gp2[n+k]
   */
  template <typename T1, std::size_t P1,
            typename T2, std::size_t P2>
  GradedPolynomial&
  operator+=(const gp_et<GradedPolynomial<T1,D,P1>, gp_div, GradedPolynomial<T2,D,P2> >& x) {
    using OuterSequence = fmmtl::GradedMonomialSequence<std::min(P,P2), DIM>;

    using first = begin_t<sequence>;
    using last  = end_t<sequence>;
    for_each<first, last>(_sum_to_mag<ORDER,MonoArray>{m_, t.m_, x.a1.m_});
    return *this;
  }

  template <typename T2, std::size_t D2, std::size_t P2>
  gp_et<GradedPolynomial,gp_mult,GradedPolynomial<T2,D2,P2> >
  operator*(const GradedPolynomial<T2,D2,P2>& gp) const {
    return {*this, gp}
  }

  template <typename T2, std::size_t D2, std::size_t P2>
  gp_et<GradedPolynomial,gp_div,GradedPolynomial<T2,D2,P2> >
  operator*(const GradedPolynomial<T2,D2,P2>& gp) const {
    return {*this, gp}
  }
#endif
};


// TODO: Specialization for DIM=1?

/** Inner product */
template <typename Array, typename T, std::size_t DIM, std::size_t ORDER>
T inner_prod(const pow_et<Array,ORDER>& x,
             const GradedPolynomial<T,DIM,ORDER>& g) {
  // Precompute the multiindex power: x^m / m!
  GradedPolynomial<T,DIM,ORDER> t(x.r, T{1});

  // TODO: Fuse with pow
  return std::inner_product(g.begin(), g.end(), t.begin(), T{});
}


template <std::size_t N, typename T, std::size_t D, std::size_t P>
const T& get(const GradedPolynomial<T,D,P>& g) {
  return get<N>(g.m_);
}

template <std::size_t N, typename T, std::size_t D, std::size_t P>
T& get(GradedPolynomial<T,D,P>& g) {
  return get<N>(g.m_);
}

}

#if 0
using namespace fmmtl;

#include <iostream>
#include <array>
#include "fmmtl/util/Clock.hpp"


int main() {
  GradedPolynomial<double,3,5> gp1, gp2;
  gp1.fill(0);
  gp2.fill(0);

  std::array<double,3> x = {2,3,5};
  constexpr multiindex<5> k;

  gp2 += pow(x,k) * 1.0;
  //gp2.fill(1);

  gp1 += pow(x,k) / gp2;

  gp1 += pow(x,k) * gp2;

  //printf("%f\n", gp2.m_[0]);


  for (auto i : gp1.m_)
    std::cout << i << std::endl;



  //for_each(fmmtl::GradedMonomialSequence<4,3>{}, printer());

  //using firsti  = typename begin<fmmtl::GradedMonomialSequence<4,3>>::type;
  //using secondi = typename next<firsti>::type;
  //using endi    = typename end<fmmtl::GradedMonomialSequence<4,3>>::type;
  //for_each(secondi{}, endi{}, printer());


  //for_each(make_index_range<0,900>{}, printer());

  //std::cout << combination(19,3) << std::endl;

  //std::cout << MultiIndex<0,0,1>::I << std::endl;
  //std::cout << monomial_index(0,0,0) << std::endl;

  /*
  std::array<double,3> x = {2,3,5};

  GradedPolynomial<double,3,5> gp(x);

  for (auto& i : gp.mono)
    std::cout << i << std::endl;
  */

  //for_each(make_index_range<0,900>{}, tester{a});

  //boost::mpl::for_each<boost::mpl::range_c<int,0,800>>( printer() );

  return 0;
}
#endif
