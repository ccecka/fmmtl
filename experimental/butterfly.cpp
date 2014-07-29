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

//#include "Fourier.kern"
//#include "ButterflyExpansion.hpp"

#include "util/Chebyshev.hpp"
#include "util/Lagrange.hpp"



#include <boost/math/special_functions/sin_pi.hpp>
#include <boost/math/special_functions/cos_pi.hpp>
#include <boost/math/constants/constants.hpp>


template <typename T>
fmmtl::complex<T> unit_polar(const T& theta) {
  using std::sin;  using std::cos;
  return fmmtl::complex<T>(cos(theta), sin(theta));
}

template <typename T>
constexpr T pow_(const T& base, const unsigned exponent) {
  return (exponent == 0) ? 1 : base * pow(base, exponent-1);
}


// TODO: Remove dependency of Direct on fmmtl::Kernel
struct FourierKernel : public fmmtl::Kernel<FourierKernel> {
  typedef double value_type;

  typedef Vec<1,value_type> source_type;
  typedef Vec<1,value_type> target_type;
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
    return 1;
    //double v = boost::math::sin_pi(t[0] - s[0]);
    //return -v * v;
    //return 2 * boost::math::constants::pi<value_type>() * (t[0] + s[0]);
  }

  value_type ampl(const target_type& t,
                  const source_type& s) const {
    (void) t; (void) s;
    return 1 + t[0] + s[0];
  }
};





/** Quickie class for iterating over the integer tensor range
 * (0,0,...,0) x (Q,Q,...,Q)    (where cardinality is DIM)
 */
template <std::size_t DIM, std::size_t Q>
struct TensorIndexGridRange {
  std::array<std::size_t,DIM> i_;

  // Current != end of the range
  template <typename T>
  bool operator!=(T&&) const {
    return i_[DIM-1] != Q;
  }
  // Increment the grid index and carry into the next dimensions
  void operator++() {
    auto it = i_.begin();
    while (++(*it) == Q && it != i_.end()-1) {
      *it = 0;
      ++it;
    }
  }
  // Return the grid index
  const std::array<std::size_t,DIM>& operator*() const {
    return i_;
  }

  TensorIndexGridRange<DIM,Q>& begin() { return *this; }
  TensorIndexGridRange<DIM,Q>& end()   { return *this; }
};





int main(int argc, char** argv) {
  int N = 10000;
  int M = 10000;
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
  typedef FourierKernel Kernel;
  Kernel K;
  double _M_ = 1;

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

  // Dimension of the tree
  const unsigned D = fmmtl::dimension<source_type>::value;
  static_assert(D == fmmtl::dimension<target_type>::value, "Dimension mismatch");

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

  // Quadrature size to use
  const unsigned Q = 10;

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
  int L_max = 2;

  std::cout << "L_max = " << L_max << std::endl;

  // Initialize all multipole/local    TODO: improve memory requirements
  for (int L = 0; L <= L_max; ++L) {
    for (source_box_type sbox : boxes(L_max - L, source_tree))
      multipole[sbox].resize(target_tree.boxes(L));
    for (target_box_type tbox : boxes(L, target_tree))
      local[tbox].resize(source_tree.boxes(L_max - L));
  }

  // Compute the level of the split
  //int L_split = L_max / 2;
  //assert(L_split > 0);

  // Begin traversal
  int L = 0;

  std::cout << "S2M" << std::endl;
  { ScopeClock timer("S2M: ");

  // For L = 0, S2M
  for (source_box_type sbox : boxes(L_max - L, source_tree)) {

    const point_type& s_center = sbox.center();
    const point_type& s_extent = sbox.extents();

    // Transform sources to reference grid
    std::vector<source_type> s_ref(p_sources[sbox.body_begin()],
                                   p_sources[sbox.body_end()]);
    for (auto&& s : s_ref) {
      s -= s_center;
      s /= s_extent;
    }

    // For each multiindex
    int i_idx = -1;
    for (auto&& i : TensorIndexGridRange<D,Q>()) {
      ++i_idx;

      // Precompute the Lagrange coefficients for the multi-index i
      std::vector<double> ls = Lagrange<Q>::coeff(i, s_ref);

      // For all the boxes in this level of the target tree
      int t_idx = -1;
      for (target_box_type tbox : boxes(L, target_tree)) {
        ++t_idx;

        const point_type& t_center = tbox.center();

        // Accumulate into M_AB_i
        complex_type& M_AB_i = multipole[sbox][t_idx][i_idx];

        auto si = p_sources[sbox.body_begin()];
        auto ci = p_charges[sbox.body_begin()];
        auto li = ls.begin();
        auto li_end = ls.end();
        for ( ; li != li_end; ++si, ++ci, ++li) {
          M_AB_i += (*li) * unit_polar(_M_ * K.phase(t_center, *si)) * (*ci);
        }
      }
    }
  }

  } // timer

#if 0
  ++L;

  std::cout << "M2M" << std::endl;

  // For all levels up to the split
  for ( ; L < L_split; ++L) {

    // For all source boxes
    for (source_box_type sbox : boxes(L_max - L, source_tree)) {

      const point_type& s_center = sbox.center();
      const point_type& s_extent = sbox.extents();

      // TOFO: Fix
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

        // Define the reference Chebyshev grid on the child box
        // TODO: encapsulate and optimize by putting into interpolant computation
        std::array<point_type,pow_(Q,D)> ref_cheb = c_cheb;
        for (auto&& r : ref_cheb) {
          r -= s_center;
          r /= s_extent;
        }

        // Accumulate
        int i_idx = -1;
        for (auto&& i : TensorIndexGridRange<D,Q>()) {
          ++i_idx;

          // Precompute Lagrange interpolants for this i
          // XXX Optimize, this can't be evaluated explicitly everytime
          std::array<double,pow_(Q,D)> ls = Lagrange<Q>::coeff(i, ref_cheb);

          // Accumulate into M_AB_t
          int t_idx = -1;
          for (target_box_type tbox : boxes(L, target_tree)) {
            ++t_idx;

            const point_type& t_center  = tbox.center();
            target_box_type pbox = tbox.parent();
            const point_type& p_center = pbox.center();
            int p_idx = pbox.index() - boxes(L-1, target_tree).begin()->index();

            // Accumulate into M_AB_i
            complex_type& M_AB_i = multipole[sbox][t_idx][i_idx];

            // For each element of i_prime
            auto li = ls.begin();
            auto li_end = ls.end();
            auto yi = c_cheb.begin();
            auto mi = multipole[cbox][p_idx].begin();
            for ( ; li != li_end; ++li, ++yi, ++mi)
              M_AB_i += (*li) * (*mi)
                  * unit_polar(_M_ * (K.phase(t_center, *yi) -
                                      K.phase(p_center, *yi)));
          }
        }
      }
    }
  }


  // For L == L_split

  // For all levels past L_split
  for ( ; L <= L_max; ++L) {

  }
#endif


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

        // TODO: Add condition for S2T...

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

  std::cout << "L2T" << std::endl;

  { ScopeClock timer("L2T: ");

  // L2T
  for (target_box_type tbox : boxes(L, target_tree)) {

    const point_type& t_center = tbox.center();
    const point_type& t_extent = tbox.extents();

    // Transform targets to reference grid
    std::vector<target_type> t_ref(p_targets[tbox.body_begin()],
                                   p_targets[tbox.body_end()]);
    for (auto&& t : t_ref) {
      t -= t_center;
      t /= t_extent;
    }

    // For each multiindex
    int i_idx = -1;
    for (auto&& i : TensorIndexGridRange<D,Q>()) {
      ++i_idx;

      // Precompute the Lagrange coefficients for the multi-index i
      std::vector<double> ls = Lagrange<Q>::coeff(i, t_ref);

      // For all the boxes in this level of the source tree
      int s_idx = -1;
      for (source_box_type sbox : boxes(L_max-L, source_tree)) {
        ++s_idx;

        const point_type& s_center = sbox.center();

        // Accumulate
        complex_type& L_AB_i = local[tbox][s_idx][i_idx];

        auto ti = p_targets[tbox.body_begin()];
        auto ri = p_results[tbox.body_begin()];
        auto li = ls.begin();
        auto li_end = ls.end();
        for ( ; li != li_end; ++ri, ++ti, ++li)
          *ri += (*li) * L_AB_i * unit_polar(_M_ * K.phase(*ti,s_center));
      }
    }
  }

  } // timer

  // Copy back hack
  auto pri = target_tree.body_permute(results.begin(),
                                      target_tree.body_begin());
  for (auto ri = p_results[target_tree.body_begin()];
       ri != p_results[target_tree.body_end()];
       ++ri, ++pri) {
    *pri += *ri;
  }

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
