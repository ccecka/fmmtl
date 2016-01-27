#pragma once
/** @file BarycentricTaylor.hpp
 * @brief Implements the Barycentric kernel with cartesian Taylor expansions.
 *
 * K(t,s) = 1 / s-t
 */

#include "fmmtl/Expansion.hpp"
#include "fmmtl/numeric/Vec.hpp"

#include "Util/Lagrange.hpp"

template <typename T>
constexpr T meta_pow(const T& base, const unsigned exponent) {
  return (exponent == 0) ? 1 : base * pow(base, exponent-1);
}



/** BBFMM -- Generic Chebyshev black-box FMM
 * @tparam P       The order of the Chebyshev polynomial expansion
 * @tparam Kernel  The base kernel to expand
 */
template <unsigned P, class Kernel>
class BBFMM
    : public fmmtl::Expansion<Kernel, BBFMM<P,Kernel> >
{
 public:
  FMMTL_IMPORT_KERNEL_TRAITS(Kernel);

  // Get the dimensions
  constexpr static int D = std::tuple_size<target_type>::value;
  static_assert(D == std::tuple_size<source_type>::value,
                "Source and Target types must have the same dimensions!");

  //! Point type
  typedef Vec<D,double> point_type;

  // Expansion type
  struct MLData {
    point_type extents;
    std::array<charge_type,meta_pow(P,D)> data;
  };

  //! Multipole expansion type
  typedef MLData multipole_type;
  //! Local expansion type
  typedef MLData local_type;

  // Base class type
  using super_t = fmmtl::Expansion<Kernel, BBFMM<P,Kernel> >;
  // Inherit base class constructors
  using super_t::super_t;

  // Initialization
  void init_multipole(multipole_type& M,
                      const point_type& extents, unsigned) const {
    M.extents = extents;
    M.data.fill(0);
  }
  void init_local(local_type& L,
                  const point_type& extents, unsigned) const {
    L.extents = extents;
    L.data.fill(0);
  }


  // Define the vector S2M
  template <class SourceIter, class ChargeIter>
  void S2M(SourceIter s_first, SourceIter s_last, ChargeIter c_first,
           const point_type& center, multipole_type& M) const {
    // Precompute Lagrange matrix with sources transformed to reference grid
    const point_type scale = 1. / M.extents;
    auto make_ref = [&](const source_type& s) {
      return (s - center) * scale;
    };
    auto LgM = LagrangeMatrix<D,P>(s_first, s_last, make_ref);

    // Multiply M_A_i += LgM_ij c_j
    prod(LgM, c_first, M.data.begin());
  }

  // Define the M2M
  void M2M(const multipole_type& Ms,
                 multipole_type& Mt,
           const point_type& translation) const {
    // D == 1 since we don't have a multidimensional Chebyshev range yet...
    // With tree structure this could be completely static too...
    static_assert(D == 1, "Only D == 1 for now...");

    // Precompute Lagrange matrix with grid transformed to reference grid
    const point_type scale = Ms.extents / Mt.extents;
    const point_type center = translation / Mt.extents;
    auto make_ref = [&](const double& c) {
      return c * scale - center;
    };
    auto LgM = LagrangeMatrix<D,P>(&Chebyshev<double,P>::x[0],
                                   &Chebyshev<double,P>::x[P],
                                   make_ref);

    // Multiply M_A_i += LgM_ij m_j
    prod(LgM, Ms.data.begin(), Mt.data.begin());
  }

  // Define the M2L
  void M2L(const multipole_type& M,
                 local_type& L,
           const point_type& translation) const {
    // Needs to be D == 1
    // Needs to be translation invariant...
    // Needs to be convertible to source_type/target_type...

    for (int i = 0; i < P; ++i) {
      target_type t = {Chebyshev<double,P>::x[i] * L.extents + translation};
      for (int j = 0; j < P; ++j) {
        source_type s = {Chebyshev<double,P>::x[j] * M.extents};
        L.data[i] += this->operator()(t, s) * M.data[j];
      }
    }
  }

  // Define the L2L
  void L2L(const local_type& Ls,
                 local_type& Lt,
           const point_type& translation) const {
    // D == 1 since we don't have a multidimensional Chebyshev range yet...
    // With tree structure this could be completely static too...
    static_assert(D == 1, "Only D == 1 for now...");

    // Precompute Lagrange matrix with grid transformed to reference grid
    const point_type scale = Lt.extents / Ls.extents;
    const point_type center = translation / Ls.extents;
    auto make_ref = [&](const double& c) {
      return c * scale + center;
    };
    auto LgM = LagrangeMatrix<D,P>(&Chebyshev<double,P>::x[0],
                                   &Chebyshev<double,P>::x[P],
                                   make_ref);

    // Multiply M_A_i += LgM_ij l_j
    prod(trans(LgM), Ls.data.begin(), Lt.data.begin());
  }

  // Define the L2T
  template <class TargetIter, class ResultIter>
  void L2T(const local_type& L, const point_type& center,
           TargetIter t_first, TargetIter t_last, ResultIter r_first) const {
    // Precompute Lagrange matrix with targets transformed to reference grid
    const point_type scale = 1. / L.extents;
    auto make_ref = [&](const target_type& t) {
      return (t - center) * scale;
    };
    auto LgM = LagrangeMatrix<D,P>(t_first, t_last, make_ref);

    // Multiply M_A_i += LgM_ij l_j
    prod(trans(LgM), L.data.begin(), r_first);
  }
};
