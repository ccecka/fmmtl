#include <iostream>

#include "fmmtl/Kernel.hpp"

#include "fmmtl/numeric/Vec.hpp"
#include "fmmtl/numeric/Complex.hpp"
#include "fmmtl/numeric/norm.hpp"

#include "fmmtl/tree/NDTree.hpp"
#include "fmmtl/tree/TreeRange.hpp"
#include "fmmtl/tree/TreeData.hpp"

#include "fmmtl/Direct.hpp"

#include "fmmtl/util/Clock.hpp"

#include "util/Chebyshev.hpp"
#include "util/Lagrange.hpp"



#include <boost/math/special_functions/sin_pi.hpp>
#include <boost/math/special_functions/cos_pi.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/iterator/transform_iterator.hpp>


//#include "ButterflyExpansion.hpp"

//! Dimension, TEMP
const std::size_t D = 2;
//! Quadrature, TEMP
const std::size_t Q = 10;

// TODO: Remove dependency of Direct on fmmtl::Kernel
struct FourierKernel : public fmmtl::Kernel<FourierKernel> {
  typedef double value_type;

  typedef Vec<D,value_type> source_type;
  typedef Vec<D,value_type> target_type;
  typedef fmmtl::complex<value_type> charge_type;
  typedef fmmtl::complex<value_type> result_type;

  typedef fmmtl::complex<value_type> kernel_value_type;

  kernel_value_type operator()(const target_type& t,
                               const source_type& s) const {
    (void) t; (void) s;
    return fmmtl::polar(ampl(t,s), phase(t,s));
    //const value_type r = 2 * (t[0] + s[0]);
    //return kernel_value_type(boost::math::cos_pi(r), boost::math::sin_pi(r));
  }

  value_type phase(const target_type& t,
                   const source_type& s) const {
    (void) t; (void) s;
    return t[0] + s[0];
    //double v = boost::math::sin_pi(t[0] - s[0]);
    //return -v * v;
    //return 2 * boost::math::constants::pi<value_type>() * (t[0] + s[0]);
  }

  value_type ampl(const target_type& t,
                  const source_type& s) const {
    (void) t; (void) s;
    return 1 + t[0] + s[0] + t[1] + s[1];
  }
};


/** Quickie class for iterating over the integer tensor range
 * (0,0,...,0) x (Q,Q,...,Q)    (where cardinality is DIM)
 */
template <std::size_t DIM, std::size_t Q>
struct TensorIndexGridRange {
  typedef std::array<unsigned,DIM> value_type;
  value_type i_ = {{}};

  // Prevent copying the value_type when iterating
  struct Reference {
    TensorIndexGridRange<DIM, Q>& t;

    // Increment the grid index and carry into the next dimensions
    void operator++() {
      for (std::size_t k = 0; ++t.i_[k] == Q && k != DIM-1; ++k)
        t.i_[k] = 0;
    }
    // Current != end of the range
    template <typename T>
    bool operator!=(T&&) const {
      return t.i_[DIM-1] != Q;
    }
    // Return the grid index
    const value_type& operator*() const {
      return t.i_;
    }
  };

  Reference begin() { return Reference{*this}; }
  Reference end()   { return Reference{*this}; }
};

template <typename T>
fmmtl::complex<T> unit_polar(const T& theta) {
  using std::sin;  using std::cos;
  return fmmtl::complex<T>(cos(theta), sin(theta));
}

template <typename T>
constexpr T pow_(const T& base, const unsigned exponent) {
  return (exponent == 0) ? 1 : base * pow(base, exponent-1);
}


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











int main(int argc, char** argv) {
  int N = 20000;
  int M = 20000;
  bool checkErrors = true;

  // Parse custom command line args
  for (int i = 1; i < argc; ++i) {
    if (strcmp(argv[i],"-N") == 0) {
      N = atoi(argv[++i]);
    } else if (strcmp(argv[i],"-M") == 0) {
      M = atoi(argv[++i]);
    } else if (strcmp(argv[i],"-nocheck") == 0) {
      checkErrors = false;
    }
  }

  // Define the kernel
  typedef ButterflyExpansion<FourierKernel, Q> Kernel;
  Kernel K;
  double _M_ = 1;  // TEMP

  std::cout << KernelTraits<Kernel>() << std::endl;

  typedef typename Kernel::source_type source_type;
  typedef typename Kernel::charge_type charge_type;
  typedef typename Kernel::target_type target_type;
  typedef typename Kernel::result_type result_type;
  typedef typename Kernel::kernel_value_type kernel_value_type;

  // Init sources
  std::vector<source_type> sources = fmmtl::random_n(N);

  // Init charges
  std::vector<charge_type> charges = fmmtl::random_n(N);

  // Init targets
  std::vector<target_type> targets = fmmtl::random_n(M);

  // Init results
  std::vector<result_type> results(M);

  // Construct two trees
  fmmtl::NDTree<D> source_tree(sources, 16);
  fmmtl::NDTree<D> target_tree(targets, 16);

  typedef typename fmmtl::NDTree<D>::box_type target_box_type;
  typedef typename fmmtl::NDTree<D>::box_type source_box_type;
  typedef typename fmmtl::NDTree<D>::body_type target_body_type;
  typedef typename fmmtl::NDTree<D>::body_type source_body_type;

  typedef typename fmmtl::NDTree<D>::point_type point_type;

  //
  // Butterfly Application
  //

  typedef fmmtl::complex<double> complex_type;

  // Associate a multipoleAB with each source box
  typedef std::vector<std::array<complex_type,pow_(Q,D)>> multipole_type;
  auto multipole = make_box_binding<multipole_type>(source_tree);

  // Associate a localAB with each target box
  typedef std::vector<std::array<complex_type,pow_(Q,D)>> local_type;
  auto local = make_box_binding<local_type>(target_tree);

  // Permute the sources and charges
  auto p_sources = make_body_binding(source_tree, sources);
  auto p_charges = make_body_binding(source_tree, charges);

  // Permute the targets and results
  auto p_targets = make_body_binding(target_tree, targets);
  auto p_results = make_body_binding(target_tree, results);

  // The maximum level of interaction
  //int L_max = std::min(source_tree.levels(), target_tree.levels()) - 1;
  int L_max = 4;
  // Compute the level of the split
  int L_split = 2;
  //assert(L_split > 0);

  std::cout << "L_max = " << L_max << std::endl;

  // Initialize all multipole/local    TODO: improve memory requirements
  for (int L = 0; L <= L_max; ++L) {
    for (source_box_type sbox : boxes(L_max - L, source_tree))
      multipole[sbox].resize(target_tree.boxes(L));
    for (target_box_type tbox : boxes(L, target_tree))
      local[tbox].resize(source_tree.boxes(L_max - L));
  }


  // Begin traversal
  int L = 0;

  std::cout << "S2M" << std::endl;
  { ScopeClock timer("S2M: ");

  // For L = 0, S2M
  for (source_box_type sbox : boxes(L_max - L, source_tree)) {
#if 1
    K.S2M(p_sources[sbox], p_charges[sbox.body_begin()],
          sbox.center(), sbox.extents(),
          // This could be better...
          boost::make_transform_iterator(target_tree.box_begin(L),
                                         [](const target_box_type& tbox) {
                                           return tbox.center();
                                         }),
          multipole[sbox]);
#else
    const point_type& s_center = sbox.center();
    const point_type s_scale = 1. / sbox.extents();

    // Precompute Lagrange interp matrix with sources transformed to ref grid
    auto make_ref = [&](const source_type& s) {
      return (s - s_center) * s_scale;
    };
    auto LgM = LagrangeMatrix<D,Q>(
        boost::make_transform_iterator(p_sources[sbox.body_begin()], make_ref),
        boost::make_transform_iterator(p_sources[sbox.body_end()],   make_ref));

    // For all the boxes in this level of the target tree
    auto tboxi = target_tree.box_begin(L);
    for (auto& M_AB : multipole[sbox]) {
      // The target center for this multipole
      const point_type& t_center = (*tboxi).center();
      ++tboxi;

      // TODO: Precompute the phase * charge

      // For each multiindex
      auto mabi = M_AB.begin();
      for (auto&& i : TensorIndexGridRange<D,Q>()) {

        // Accumulate into M_AB_i
        complex_type& M_AB_i = *mabi;
        ++mabi;

        auto ci = p_charges[sbox.body_begin()];
        std::size_t j = 0;   // XXX, Abstract as mat-vec
        for (auto&& s : p_sources[sbox]) {
          M_AB_i += unit_polar(_M_ * K.phase(t_center, s)) * LgM(i,j) * (*ci);
          ++ci; ++j;
        }
      }
    }
#endif
  }

  } // timer

  ++L;

  std::cout << "M2M" << std::endl;
  { ScopeClock timer("M2M: ");

  // For all levels up to the split
  for ( ; L <= L_split; ++L) {

    // For all source boxes
    for (source_box_type sbox : boxes(L_max - L, source_tree)) {

      const point_type& s_center = sbox.center();
      const point_type  s_scale  = 1. / sbox.extents();

      // TODO: Fix
      assert(!sbox.is_leaf());

      // For all the children of this source box
      for (source_box_type cbox : children(sbox)) {

        const point_type& c_center = cbox.center();
        const point_type& c_extent = cbox.extents();

        // Define the child box chebyshev grid
        std::array<point_type,pow_(Q,D)> c_cheb;
        auto cbi = c_cheb.begin();
        for (auto&& i : TensorIndexGridRange<D,Q>()) {
          for (unsigned d = 0; d < D; ++d)
            (*cbi)[d] = c_center[d]
                + Chebyshev<double,Q>::x[i[d]] * c_extent[d];
          ++cbi;
        }

        // Precompute Lagrange interp matrix with points transformed to ref grid
        // XXX: Avoid computing at all for quad/octrees
        auto make_ref = [&](const point_type& c) {
          return (c - s_center) * s_scale;
        };
        auto LgM = LagrangeMatrix<D,Q>(
            boost::make_transform_iterator(c_cheb.begin(), make_ref),
            boost::make_transform_iterator(c_cheb.end(),   make_ref));

        // Accumulate into M_AB_t
        int t_idx = -1;
        for (target_box_type tbox : boxes(L, target_tree)) {
          ++t_idx;

          const point_type& t_center  = tbox.center();
          target_box_type pbox = tbox.parent();
          const point_type& p_center = pbox.center();
          int p_idx = pbox.index() - boxes(pbox.level(), target_tree).begin()->index();

          // Accumulate
          int i_idx = -1;
          for (auto&& i : TensorIndexGridRange<D,Q>()) {
            ++i_idx;

            // Accumulate into M_AB_i
            complex_type& M_AB_i   = multipole[sbox][t_idx][i_idx];

            // For each element of i_prime
            auto yi = c_cheb.begin();
            std::size_t j = 0;  // Lift the matvec
            for (auto&& M_ApBc_ip : multipole[cbox][p_idx]) {
              M_AB_i += M_ApBc_ip  * LgM(i,j)
                  * unit_polar(_M_ * (K.phase(t_center, *yi) -
                                      K.phase(p_center, *yi)));
              ++yi; ++j;
            }
          }
        }
      }
    }
  }

  }

  // M2L on the last level that was M2Med
  --L;

  std::cout << "M2L" << std::endl;
  { ScopeClock timer("M2L: ");

  // M2L and check the result
  int s_idx = -1;
  for (source_box_type sbox : boxes(L_max - L, source_tree)) {
    ++s_idx;

    const point_type& s_center = sbox.center();
    const point_type& s_extent = sbox.extents();

    // Define the source box Chebyshev grid
    std::array<point_type,pow_(Q,D)> s_cheb;
    {
    auto si = s_cheb.begin();
    for (auto&& i : TensorIndexGridRange<D,Q>()) {
      for (unsigned d = 0; d < D; ++d)
        (*si)[d] = s_center[d]
            + Chebyshev<double,Q>::x[i[d]] * s_extent[d];
      ++si;
    }
    }

    int t_idx = -1;
    for (target_box_type tbox : boxes(L, target_tree)) {
      ++t_idx;

      const point_type& t_center = tbox.center();
      const point_type& t_extent = tbox.extents();

      // Define the target box Chebyshev grid
      std::array<point_type,pow_(Q,D)> t_cheb;
      {
      auto ti = t_cheb.begin();
      for (auto&& i : TensorIndexGridRange<D,Q>()) {
        for (unsigned d = 0; d < D; ++d)
          (*ti)[d] = t_center[d]
              + Chebyshev<double,Q>::x[i[d]] * t_extent[d];
        ++ti;
      }
      }

      auto li = local[tbox][s_idx].begin();
      auto ti = t_cheb.begin();
      auto ti_end = t_cheb.end();
      for ( ; ti != ti_end; ++ti, ++li) {

        // Grab the multipole
        auto mi = multipole[sbox][t_idx].begin();
        auto si = s_cheb.begin();
        auto si_end = s_cheb.end();
        for ( ; si != si_end; ++si, ++mi) {
          *li += K.ampl(*ti,*si) * (*mi)
              * unit_polar(_M_ * (K.phase(*ti,*si) - K.phase(t_center,*si)));
        }
        *li *= unit_polar(-_M_ * K.phase(*ti, s_center));
      }
    }
  }

  } // timer


  std::cout << "L2L" << std::endl;
  { ScopeClock timer("L2L: ");

  // For all levels up to the max
  for ( ; L <= L_max; ++L) {

    // For all target boxes
    for (target_box_type tbox : boxes(L, target_tree)) {

      const point_type& t_center = tbox.center();
      const point_type& t_extent = tbox.extents();

      // Define the target box Chebyshev grid
      std::array<point_type,pow_(Q,D)> t_cheb;
      auto cbi = t_cheb.begin();
      for (auto&& i : TensorIndexGridRange<D,Q>()) {
        for (unsigned d = 0; d < D; ++d)
          (*cbi)[d] = t_center[d]
              + Chebyshev<double,Q>::x[i[d]] * t_extent[d];
        ++cbi;
      }

      // Get the target box's parent
      target_box_type pbox = tbox.parent();

      const point_type& p_center = pbox.center();
      const point_type  p_scale  = 1. / pbox.extents();

      // Precompute Lagrange interp matrix with points transformed to ref grid
      // XXX: Avoid computing at all for quad/octrees
      auto make_ref = [&](const point_type& c) {
        return (c - p_center) * p_scale;
      };
      auto LgM = LagrangeMatrix<D,Q>(
          boost::make_transform_iterator(t_cheb.begin(), make_ref),
          boost::make_transform_iterator(t_cheb.end(),   make_ref));

      // For all the child boxes
      int c_idx = -1;
      for (source_box_type cbox : boxes(L_max - L + 1, source_tree)) {
        ++c_idx;

        const point_type& c_center = cbox.center();

        // Get the parent
        source_box_type sbox = cbox.parent();
        int s_idx = sbox.index() - boxes(L_max - L, source_tree).begin()->index();
        const point_type& s_center = sbox.center();

        // Accumulate
        int i_idx = -1;
        for (auto&& i : TensorIndexGridRange<D,Q>()) {
          ++i_idx;

          complex_type& L_ApBc_ip = local[pbox][c_idx][i_idx];

          // For each element of i_prime
          auto xi = t_cheb.begin();
          std::size_t j = 0;  // Lift the matvec
          for (auto&& L_AB_i : local[tbox][s_idx]) {
            L_AB_i += L_ApBc_ip * LgM(i,j)
                * unit_polar(_M_ * (K.phase(*xi, c_center) -
                                    K.phase(*xi, s_center)));
            ++xi;
            ++j;
          }
        }
      }
    }
  }

  } // timer

  --L;


  std::cout << "L2T" << std::endl;
  { ScopeClock timer("L2T: ");

  // L2T
  for (target_box_type tbox : boxes(L, target_tree)) {

    const point_type& t_center = tbox.center();
    const point_type  t_scale  = 1. / tbox.extents();

    // Precompute Lagrange interp matrix with targets transformed to ref grid
    auto make_ref = [&](const target_type& t) {
      return (t - t_center) * t_scale;
    };
    auto LgM = LagrangeMatrix<D,Q>(
        boost::make_transform_iterator(p_targets[tbox.body_begin()], make_ref),
        boost::make_transform_iterator(p_targets[tbox.body_end()],   make_ref));

    // For all the boxes in this level of the source tree
    auto sbi = source_tree.box_begin(L_max - L);
    for (const auto& L_AB : local[tbox]) {
      // The source center of L_AB
      const point_type& s_center = (*sbi).center();
      ++sbi;

      // For each multiindex
      auto labi = L_AB.begin();
      for (auto&& i : TensorIndexGridRange<D,Q>()) {

        // Compute from L_AB_i
        const complex_type& L_AB_i = *labi;
        ++labi;

        auto ri = p_results[tbox.body_begin()];
        std::size_t j = 0;     // XXX, Abstract as mat-vec
        for (auto&& t : p_targets[tbox]) {
          *ri += LgM(i,j) * L_AB_i * unit_polar(_M_ * K.phase(t, s_center));
          ++j; ++ri;
        }
      }
    }
  }

  } // timer

  // Copy back permuted
  auto pri = target_tree.body_permute(results.begin(), target_tree.body_begin());
  for (auto ri : p_results) {
    *pri = ri;
    ++pri;
  }



#if 0
  // For all levels up to the split
  for (int L = 0; L <= L_max; ++L) {

    // For all boxes in the opposing level of the source tree
    int s_idx = -1;
    for (source_box_type sbox : boxes(L_max - L, source_tree)) {
      ++s_idx;

      // For all the boxes in this level of the target tree
      int t_idx = -1;
      for (target_box_type tbox : boxes(L, target_tree)) {
        ++t_idx;

        // TODO: Add condition for S2T?

        if (L == 0 || (L <= L_split && sbox.is_leaf())) {
          // sbox is a leaf with a multipole
          // S2M
        } else if (L <= L_split) {
          // sbox has a multipole and children
          // M2M
        } else if (sbox.is_leaf()) {
          // sbox is a leaf without a multipole
          // S2L
        }

        if (L == L_split) {
          // M2L
        }

        if (L == L_max || (L >= L_split && tbox.is_leaf())) {
          // tbox is a leaf with a local
          // L2T
        } else if (L >= L_split) {
          // tbox has a local and children
          // L2L
        } else if (tbox.is_leaf()) {
          // tbox is a leaf without a local
          // M2T
        }
      }
    }
  }
#endif





  // Check the result
  if (checkErrors) {
    std::cout << "Computing direct matvec..." << std::endl;

    // Compute the result with a direct matrix-vector multiplication
    std::vector<result_type> exact(M);
    { ScopeClock timer("Direct: ");
    Direct::matvec(K, sources, charges, targets, exact);
    }

    double tot_error_sq = 0;
    double tot_norm_sq = 0;
    double tot_ind_rel_err = 0;
    double max_ind_rel_err = 0;
    for (unsigned k = 0; k < results.size(); ++k) {
      //std::cout << results[k] << "\t" << exact[k] << std::endl;

      // Individual relative error
      double rel_error = norm_2(exact[k] - results[k]) / norm_2(exact[k]);
      tot_ind_rel_err += rel_error;
      // Maximum relative error
      max_ind_rel_err  = std::max(max_ind_rel_err, rel_error);

      // Total relative error
      tot_error_sq += norm_2_sq(exact[k] - results[k]);
      tot_norm_sq  += norm_2_sq(exact[k]);
    }
    double tot_rel_err = sqrt(tot_error_sq/tot_norm_sq);
    std::cout << "Vector  relative error: " << tot_rel_err << std::endl;

    double ave_rel_err = tot_ind_rel_err / results.size();
    std::cout << "Average relative error: " << ave_rel_err << std::endl;

    std::cout << "Maximum relative error: " << max_ind_rel_err << std::endl;
  }

  // TODO: Interpolative decompositions

  // TODO: Range based iterations for tree boxes -- simplification
  // TODO: Generalizations on Context so I can use it in this mock up.
}
