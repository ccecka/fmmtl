#include <array>

#include "fmmtl/Expansion.hpp"

#include "util/Chebyshev.hpp"
#include "util/Lagrange.hpp"

/**
 * @tparam Kernel The kernel to create a butterfly expansion for.
 * Special requirements:
 *   std::is_same<typename Kernel::kernel_value_type, fmmtl::complex<double>>
 *   Kernel.phase(Vec<D,double>&, Vec<D,double>&)
 *   Kernel.ampl(Vec<D,double>&, Vec<D,double>&)
 * @tparam Q    The number of quadrature points in each dimension
 */
// TODO: Replace Expansion hierarchy with Expansion triplet for generality
template <class Kernel, std::size_t Q>
struct ButterflyExpansion : public Kernel {
  FMMTL_IMPORT_KERNEL_TRAITS(Kernel);

  static constexpr std::size_t D = fmmtl::dimension<source_type>::value;
  static_assert(D == fmmtl::dimension<target_type>::value,
                "Dimension mismatch not supported");

  typedef Vec<D,double> point_type;

  typedef fmmtl::complex<double> complex_type;

  typedef std::array<complex_type,pow_(Q,D)> multipole_type;
  typedef std::array<complex_type,pow_(Q,D)> local_type;

  double _M_ = 1;
  ButterflyExpansion(double M = 1) : _M_(M) {}


  template <typename SourceRange, typename ChargeIter,
            typename TCenterIter, typename MultipoleRange>
  void S2M(SourceRange&& s_range, ChargeIter&& c_begin,
           const point_type& s_center, const point_type& s_extent,
           TCenterIter&& tc_begin, MultipoleRange&& m_range) {

    const point_type s_scale = 1. / s_extent;

    // Precompute Lagrange interp matrix with sources transformed to ref grid
    auto make_ref = [&](const source_type& s) {
      return (s - s_center) * s_scale;
    };
    auto LgM = LagrangeMatrix<D,Q>(
        boost::make_transform_iterator(s_range.begin(), make_ref),
        boost::make_transform_iterator(s_range.end(),   make_ref));

    // For all the boxes in this level of the target tree
    for (auto& M_AB : m_range) {
      // The target center for this multipole
      const point_type& t_center = *tc_begin;
      ++tc_begin;

      // TODO: Precompute the phase * charge

      // For each multiindex
      // TODO: Encapsulate in matvec
      auto mabi = M_AB.begin();
      for (auto&& i : TensorIndexGridRange<D,Q>()) {

        // Accumulate into M_AB_i
        complex_type& M_AB_i = *mabi;
        ++mabi;

        auto ci = c_begin;
        std::size_t j = 0;   // XXX, Abstract as mat-vec
        for (auto&& s : s_range) {
          M_AB_i += unit_polar(_M_ * this->phase(t_center, s)) * LgM(i,j) * (*ci);
          ++ci; ++j;
        }
      }
    }
  }

  void M2M(const multipole_type& source,
           const point_type& s_center, const point_type& s_extent,
                 multipole_type& target,
           const point_type& t_center, const point_type& t_extent) const {

  }

  void M2L(const multipole_type& source,
                 local_type& target,
           const point_type& translation) const {
    (void) source;
    (void) target;
    (void) translation;
  }

  void L2L(const local_type& source,
                 local_type& target,
           const point_type& translation) const {
    (void) source;
    (void) target;
    (void) translation;
  }

  template <typename TargetIter, typename ResultIter>
  void M2T(const multipole_type& M, const point_type& center,
           TargetIter t_first, TargetIter t_last,
           ResultIter r_first) const {
    (void) M;
    (void) center;
    (void) t_first;
    (void) t_last;
    (void) r_first;
  }

  template <typename TargetIter, typename ResultIter>
  void L2T(const local_type& L, const point_type& center,
           TargetIter t_first, TargetIter t_last,
           ResultIter r_first) const {
    (void) L;
    (void) center;
    (void) t_first;
    (void) t_last;
    (void) r_first;
  }
};
