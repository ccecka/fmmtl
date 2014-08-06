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
#include "util/TensorIndexGridRange.hpp"


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

template <typename T>
fmmtl::complex<T> unit_polar(const T& theta) {
  using std::sin;  using std::cos;
  return fmmtl::complex<T>(cos(theta), sin(theta));
}

template <typename T>
constexpr T pow_(const T& base, const unsigned exponent) {
  return (exponent == 0) ? 1 : base * pow(base, exponent-1);
}




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
  //typedef ButterflyExpansion<FourierKernel, Q> Kernel;
  typedef FourierKernel Kernel;
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

  std::cout << "Tree Construction" << std::endl;
  Clock timer;

  // Construct two trees
  fmmtl::NDTree<D> source_tree(sources, 16);
  fmmtl::NDTree<D> target_tree(targets, 16);

  std::cout << "Tree Construction: " << timer.seconds() << std::endl;

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

  // HACKY way to define S2M, M2M, M2L, L2L, L2T
#include "butterfly_operators.hpp"


#if 0
  for (int L = 0; L <= L_max; ++L) {

    if (L == 0) {
      for (source_box_type sbox : boxes(L_max - L, source_tree))
        S2M(L, sbox);
    } else if (L <= L_split) {
      for (source_box_type sbox : boxes(L_max - L, source_tree)) {
        if (sbox.is_leaf())
          S2M(L, sbox);
        else
          M2M(L, sbox);
      }
    } else {
      for (source_box_type sbox : boxes(L_max - L, source_tree)) {
        if (sbox.is_leaf()) {
          std::cerr << "Unimplemented" << std::endl;
          exit(0);
          // S2M(L, sbox);
          // M2L(L, sbox);  // need sbox parameter
          // Replace with S2L
        }
      }
    }

    if (L == L_split)
      M2L(L);

    if (L == L_max) {
      for (target_box_type tbox : boxes(L, target_tree)) {
        L2T(L, tbox);
      }
    } else if (L >= L_split) {
      for (target_box_type tbox : boxes(L, target_tree)) {
        if (tbox.is_leaf()) {
          L2T(L, tbox);
        } else {
          L2L(L, tbox);
        }
      }
    } else {
      for (target_box_type tbox : boxes(L, target_tree)) {
        if (tbox.is_leaf()) {
          std::cerr << "Unimplemented M2T" << std::endl;
          exit(0);
          // Replace with M2T
        }
      }
    }
  }
#endif

#if 1
  // Begin traversal
  int L = 0;

  std::cout << "S2M" << std::endl;
  { ScopeClock timer("S2M: ");

  // For L = 0, S2M
  for (source_box_type sbox : boxes(L_max - L, source_tree)) {
    S2M(L, sbox);
  }

  } // timer


  std::cout << "M2M" << std::endl;
  { ScopeClock timer("M2M: ");

  // For all levels up to the split
  for (++L; L <= L_split; ++L) {

    // For all child boxes
    for (source_box_type sbox : boxes(L_max - L, source_tree)) {
      M2M(L, sbox);
    }

  }

  }

  // M2L on the last level that was M2Med
  --L;

  std::cout << "M2L" << std::endl;
  { ScopeClock timer("M2L: ");

  M2L(L);

  } // timer

  //++L;

  std::cout << "L2L" << std::endl;
  { ScopeClock timer("L2L: ");

  // For all levels up to the max
  for ( ; L < L_max; ++L) {

    // For all target boxes
    for (target_box_type tbox : boxes(L, target_tree)) {
      L2L(L, tbox);
    }
  }

  } // timer


  std::cout << "L2T" << std::endl;
  { ScopeClock timer("L2T: ");

  // L2T
  for (target_box_type tbox : boxes(L, target_tree)) {
    L2T(L, tbox);
  }

  } // timer
#endif


  // Copy back permuted
  auto pri = target_tree.body_permute(results.begin(), target_tree.body_begin());
  for (auto ri : p_results) {
    *pri = ri;
    ++pri;
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
