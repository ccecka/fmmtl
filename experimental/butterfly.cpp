#include <iostream>
#include <complex>

#include <boost/math/special_functions/sin_pi.hpp>
#include <boost/math/special_functions/cos_pi.hpp>
using namepsace boost::math;

#include "fmmtl/numeric/Vec.hpp"
#include "fmmtl/numeric/norm.hpp"

#include "fmmtl/tree/NDTree.hpp"

struct FourierKernel {
  typedef double value_type;

  typedef Vec<1,value_type> source_type;
  typedef Vec<1,value_type> target_type;
  typedef std::complex<value_type> charge_type;
  typedef std::complex<value_type> result_type;

  typedef std::complex<value_type> kernel_value_type;

  kernel_value_type operator()(const target_type& t, const source_type& s) {
    const value_type r = 2*inner_prod(t,s);
    return kernel_value_type(cos_pi(r), sin_pi(r));
  }
};


int main() {
  int N = 10000;
  int M = 10000;

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

  typedef typename Kernel::source_type source_type;
  typedef typename Kernel::charge_type charge_type;
  typedef typename Kernel::target_type target_type;
  typedef typename Kernel::result_type result_type;

  // Init sources
  std::vector<source_type> sources(N);
  for (source_type& s : sources)
    s = fmmtl::random<source_type>::get();

  // Init charges
  std::vector<charge_type> charges(N);
  for (charge_type& c : charges)
    c = fmmtl::random<charge_type>::get();

  // Init targets
  std::vector<target_type> targets(M);
  for (target_type& t : targets)
    t = fmmtl::random<target_type>::get();

  // Dimension of the tree
  const unsigned D = fmmtl::dimension<source_type>::value;
  static_assert(D == fmmtl::dimension<target_type>::value);

  // Construct two trees
  fmmtl::NDTree<D> source_tree(sources);
  fmmtl::NDTree<D> target_tree(targets);

  typedef typename fmmtl::NDTree<D>::box_type target_box_type;
  typedef typename fmmtl::NDTree<D>::box_type source_box_type;

  //
  // Butterfly precomputation
  //

  // For all the levels of the target_tree
  int max_level = std::min(source_tree.levels(), target_tree.levels());
  for (int level = 0; level < max_level; ++level) {

    // For all the boxes in this level of the target tree
    auto tb_end = target_tree.box_begin(level);
    for (auto tbi = target_tree.box_begin(level); tbi != tb_end; ++tbi) {
      target_box_type tbox = *tbi;

      // For all boxes in the opposing level of the source tree
      int slevel = max_level - level - 1;
      assert(slevel >= 0);
      auto sb_end = source_tree.box_begin(slevel);
      for (auto sbi = source_tree.box_begin(slevel); sbi != sb_end; ++sbi) {
        source_box_type sbox = *sbi;

        // Interpolative decomposition of tbox x sbox
        // or propogation of previous interpolative decompositions
      }
    }
  }

  //
  // Butterfly Application
  //

  // For all the levels of the target_tree
  int max_level = std::min(source_tree.levels(), target_tree.levels());
  for (int level = 0; level < max_level; ++level) {

    // For all the boxes in this level of the target tree
    auto tb_end = target_tree.box_begin(level);
    for (auto tbi = target_tree.box_begin(level); tbi != tb_end; ++tbi) {
      target_box_type tbox = *tbi;

      // For all boxes in the opposing level of the source tree
      int slevel = max_level - level - 1;
      assert(slevel >= 0);
      auto sb_end = source_tree.box_begin(slevel);
      for (auto sbi = source_tree.box_begin(slevel); sbi != sb_end; ++sbi) {
        source_box_type sbox = *sbi;

        // Application of precomputed interpolative decompositions
      }
    }
  }


  // TODO: Interpolative decompositions

  // TODO: Range based iterations for tree boxes -- simplification
  // TODO: Generalizations on Context so I can use it in this mock up.
}
