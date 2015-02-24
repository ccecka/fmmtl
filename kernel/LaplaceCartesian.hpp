#pragma once
/** @file LaplaceCartesian.hpp
 * @brief Implements the Laplace kernel with Taylor expansions.
 *
 * K(t,s) = 1 / |s-t|        // Laplace potential
 * K(t,s) = (s-t) / |s-t|^3  // Laplace force
 */

#include "fmmtl/Expansion.hpp"
#include "fmmtl/numeric/Vec.hpp"

#include "Laplace.kern"
#include "Util/GradedPolynomial.hpp"
#include "fmmtl/meta/for_each.hpp"

/** Recurrence for computing the gradient of the 3D Laplace potential
 *
 * Let S_n = Del^n 1/||r||
 * where n = (n_0,n_1,n_2) and r = (r_0,r_1,r_2), then
 *
 * S_n = 1/||r||^2 * [(1-2|n|)/|n| Sum_i n_i r_i S_{n-e_i} +
 *                    (1- |n|)/|n| Sum_i n_i (n_i-1) S_{n-2e_i}]
 * where |n|    = n_0 + n_1 + n_2
 *      ||r||^2 = r_0^2 + r_1^2 + r_2^2
 *        e_0   = (1,0,0), e_1   = (0,1,0), e_2   = (0,0,1)
 */
template <typename Array>
struct LaplaceDerivativeRecurrence {
  template <typename N>
  void operator()(N) {
    using e0 = fmmtl::MultiIndex<1,0,0>;
    using e1 = fmmtl::MultiIndex<0,1,0>;
    using e2 = fmmtl::MultiIndex<0,0,1>;
    using std::get;

    constexpr double c1 = (1.0-2*N::R) / N::R;
    constexpr double c2 = (1.0-  N::R) / N::R;
    constexpr std::size_t n0 = get<0>(N{});
    constexpr std::size_t n1 = get<1>(N{});
    constexpr std::size_t n2 = get<2>(N{});

    double sn = 0;
    if (n0 > 0)   // static evaluation of if-statements
      sn += (c1*n0) * get<0>(r) * get<decltype(N{}-e0{})::I>(s);
    if (n1 > 0)
      sn += (c1*n1) * get<1>(r) * get<decltype(N{}-e1{})::I>(s);
    if (n2 > 0)
      sn += (c1*n2) * get<2>(r) * get<decltype(N{}-e2{})::I>(s);

    if (n0 > 1)
      sn += (c2*n0*(n0-1)) * get<decltype(N{}-e0{}-e0{})::I>(s);
    if (n1 > 1)
      sn += (c2*n1*(n1-1)) * get<decltype(N{}-e1{}-e1{})::I>(s);
    if (n2 > 1)
      sn += (c2*n2*(n2-1)) * get<decltype(N{}-e2{}-e2{})::I>(s);

    get<N::I>(s) = inv_r_sq * sn;
  }
  const Vec<3,double>& r;   // (r0,r1,r2)
  double inv_r_sq;          // 1.0 / |r|^2
  Array& s;
};

/** Transform the Taylor expansion into the result. In this case, the result is
 * the zeroth and first derivatives of the Laplace kernel.
 */
template <std::size_t P>
struct LaplaceL2T {
  template <typename K>
  void operator()(K) {
    using e0 = fmmtl::MultiIndex<1,0,0>;
    using e1 = fmmtl::MultiIndex<0,1,0>;
    using e2 = fmmtl::MultiIndex<0,0,1>;
    using std::get;

    result[0] += get<K::I>(r) * get<K::I>(L);
    if (get<0>(K{}) > 0)
      result[1] += get<decltype(K{}-e0{})::I>(r) * get<K::I>(L);
    if (get<1>(K{}) > 0)
      result[2] += get<decltype(K{}-e1{})::I>(r) * get<K::I>(L);
    if (get<2>(K{}) > 0)
      result[3] += get<decltype(K{}-e2{})::I>(r) * get<K::I>(L);
  }
  Vec<4,double>& result;
  const fmmtl::GradedPolynomial<double,3,P>& r;
  const fmmtl::GradedPolynomial<double,3,P>& L;
};



template <unsigned P>
class LaplaceCartesian
    : public fmmtl::Expansion<LaplaceKernel, LaplaceCartesian<P>>
{
 public:
  FMMTL_IMPORT_KERNEL_TRAITS(LaplaceKernel);

  //! Point type
  typedef Vec<3,double> point_type;

  //! Multipole expansion type
  typedef fmmtl::GradedPolynomial<double,3,P> multipole_type;
  //! Local expansion type
  typedef fmmtl::GradedPolynomial<double,3,P> local_type;

  //! A multiindex ranging from 0 to P
  static constexpr fmmtl::multiindex<P> k{};

  void init_multipole(multipole_type& M, const point_type&, unsigned) const {
    M.fill(0);
  }
  void init_local(local_type& L, const point_type&, unsigned) const {
    L.fill(0);
  }


  /** Kernel S2M operation */
  void S2M(const source_type& source, const charge_type& charge,
           const point_type& center, multipole_type& M) const {
    M += pow(center-source, k) * charge;
  }

  /** Kernel M2M operator */
  void M2M(const multipole_type& Ms,
                 multipole_type& Mt,
           const point_type& translation) const {
    Mt += pow(translation, k) * Ms;
  }

  /** Kernel M2L operation */
  void M2L(const multipole_type& M,
                 local_type& L,
           const point_type& r) const {
    // Storage for the multi-index derivative of 1/|r|
    constexpr std::size_t DelP = P;
    //constexpr std::size_t DelP = 2*P;  // Can help error, but very slow.
    using DelPoly = fmmtl::GradedPolynomial<double,3,DelP>;
    using Recurrence = LaplaceDerivativeRecurrence<DelPoly>;
    DelPoly del_f;

    // Initialize the monopole
    double invR2 = 1./norm_2_sq(r);
    del_f[0] = std::sqrt(invR2);

    // Use the recurrence to compute the derivatives
    using SequenceDel = typename DelPoly::sequence;
    using second2p    = fmmtl::next_t<fmmtl::begin_t<SequenceDel>>;
    using last2p      = fmmtl::end_t<SequenceDel>;
    fmmtl::for_each<second2p, last2p>(Recurrence{r, invR2, del_f});

    // Compute L[n] += sum_k del_f[n+k] * M[k]
    using sequence = fmmtl::GradedMonomialSequence<P,3>;
    using first    = fmmtl::begin_t<sequence>;
    using last     = fmmtl::end_t<sequence>;
    fmmtl::for_each<first, last>(fmmtl::sum_to_npk_mag<P>(L,M,del_f));
    // For use with DelP = 2*P:
    //fmmtl::for_each<first, last>(fmmtl::sum_to_k_mag<P>(L,M,del_f));
  }

  /** Kernel L2L operator */
  void L2L(const local_type& Ls,
                 local_type& Lt,
           const point_type& translation) const {
    Lt += pow(translation, k) / Ls;
  }

  /** Kernel L2T operator */
  void L2T(const local_type& L, const point_type& center,
           const target_type& target, result_type& result) const {
    // Precompute the multiindex power: r^m / m!
    fmmtl::GradedPolynomial<double,3,P> t(target-center, 1);

    using sequence = fmmtl::GradedMonomialSequence<P,3>;
    using first    = fmmtl::begin_t<sequence>;
    using last     = fmmtl::end_t<sequence>;

    fmmtl::for_each<first, last>(LaplaceL2T<P>{result,t,L});
  }
};
