#include <iostream>
#include <complex>

#include <boost/math/special_functions/sin_pi.hpp>
#include <boost/math/special_functions/cos_pi.hpp>
using namespace boost::math;

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
using namespace boost::numeric::ublas;

#include "fmmtl/Kernel.hpp"

#include "fmmtl/numeric/Vec.hpp"
#include "fmmtl/numeric/Complex.hpp"
#include "fmmtl/numeric/norm.hpp"

#include "fmmtl/tree/NDTree.hpp"
#include "fmmtl/tree/TreeRange.hpp"
#include "fmmtl/tree/TreeData.hpp"

#include "fmmtl/executor/Traversal.hpp"

#include "fmmtl/Direct.hpp"


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
    const value_type r = 2*inner_prod(t,s);
    return kernel_value_type(cos_pi(r), sin_pi(r));
  }
};


int main(int argc, char** argv) {
  int N = 1000;
  int M = 1000;
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

  std::cout << KernelTraits<Kernel>() << std::endl;

  typedef typename Kernel::source_type source_type;
  typedef typename Kernel::charge_type charge_type;
  typedef typename Kernel::target_type target_type;
  typedef typename Kernel::result_type result_type;
  typedef typename Kernel::kernel_value_type kernel_value_type;

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

  // Init results
  std::vector<result_type> result(M);

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

  // Operators associated with the source and target boxes
  auto p_op = make_box_binding<std::vector<matrix<double>>>(source_tree);
  auto b_op = make_box_binding<std::vector<matrix<kernel_value_type>>>(target_tree);

  //
  // Butterfly precomputation
  //

  // For all the levels of the target_tree
  int max_L = std::min(source_tree.levels(), target_tree.levels()) - 1;
  //int max_L = 3;

  // if maxL == L, apply direct evaluation

  std::cout << "max_L = " << max_L << std::endl;
  for (int L = 0; L <= max_L; ++L) {

    std::cout << "Level " << L << std::endl;

    // For all boxes in the opposing level of the source tree
    for (source_box_type sbox : boxes(max_L - L, source_tree)) {

      p_op[sbox].reserve(target_tree.boxes(L));

      // For all the boxes in this level of the target tree
      for (target_box_type tbox : boxes(L, target_tree)) {

        std::cout << "Precomputing " << L << "  " << tbox.index() << "  " << sbox.index() << std::endl;

        // Interpolative decomposition of tbox x sbox
        // or propogation of previous interpolative decompositions

        if (L == 0 || sbox.is_leaf()) {
          // Form interpolative decomposition from sources/targets (dummy)
          p_op[sbox].push_back(identity_matrix<double>(sbox.num_bodies(),
                                                       sbox.num_bodies()));
          //t_op[tbox.index()] = B-matrix... ?
        } else {
          // Get the "B"-matrix of the target's parent and partition according to child partitions
          // Form the interpolative decomposition of each partition (dummy)
          p_op[sbox].push_back(identity_matrix<double>(sbox.num_bodies(),
                                                       sbox.num_bodies()));

          // Since we've only done dummy P-matrices, this will be the full evaluated rows
        }

        if (L == max_L || tbox.is_leaf()) {
          b_op[tbox].reserve(source_tree.boxes(max_L - L));

          // Construct the elements of the kernel matrix
          // Could use a ublas kernel_matrix class here...

          b_op[tbox].push_back(matrix<kernel_value_type>(tbox.num_bodies(),
                                                         sbox.num_bodies()));
          auto& k = b_op[tbox].back();

          std::cout << "Matrix " << k.size1() << "," << k.size2() << std::endl;
          int tb_idx = -1;
          for (target_body_type tb : bodies(tbox)) {
            ++tb_idx;
            target_type& t = targets[tb.number()];

            int sb_idx = -1;
            for (source_body_type sb : bodies(sbox)) {
              ++sb_idx;
              source_type& s = sources[sb.number()];
              k(tb_idx, sb_idx) = K(t, s);
            }
          }
        }
      }
    }
  }

  //
  // Butterfly Application
  //

  typedef std::vector<vector<charge_type>> multipole_type;
  auto multipole = make_box_binding<multipole_type>(source_tree);

  // For all the levels of the target_tree
  for (int L = 0; L <= max_L; ++L) {

    // For all boxes in the opposing level of the source tree
    int s_idx = -1;
    for (source_box_type sbox : boxes(max_L - L, source_tree)) {
      ++s_idx;

      multipole[sbox].resize(target_tree.boxes(L));

      // For all the boxes in this level of the target tree
      int t_idx = -1;
      for (target_box_type tbox : boxes(L, target_tree)) {
        ++t_idx;

        std::cout << "Computing " << L << "  " << tbox.index() << "  " << sbox.index() << std::endl;

        // Application of precomputed interpolative decompositions
        if (L == 0 || sbox.is_leaf()) {
          // Create the "multipole" from the sources (S2M)
          // Construct the local source vector explictly...
          vector<charge_type> c(sbox.num_bodies());
          int sb_idx = -1;
          for (source_body_type sb : bodies(sbox)) {
            ++sb_idx;
            c[sb_idx] = charges[sb.number()];
          }
          // Compute the product (dummy)
          multipole[sbox][t_idx] = prod(p_op[sbox][t_idx], c);
        } else {
          // Concatenate the "multipoles" from the children and apply op (M2M)
          vector<charge_type> c(sbox.num_bodies());
          auto ci = c.begin();

          int tp_idx = tbox.parent().index() - target_tree.box_begin(L-1)->index();
          for (source_box_type cbox : children(sbox)) {
            auto& M = multipole[cbox][tp_idx];
            ci = std::copy(M.begin(), M.end(), ci);
          }
          // And take the product
          multipole[sbox][t_idx] = prod(p_op[sbox][t_idx], c);
        }

        if (L == max_L || tbox.is_leaf()) {
          std::cout << "Local Product " << tbox.index() << " - "
                    << b_op[tbox][s_idx].size1() << "x" << b_op[tbox][s_idx].size2()
                    << " . " << multipole[sbox][t_idx].size() << std::endl;
          // Take the product of the ID and the "multipole" (M2L)
          vector<result_type> c = prod(b_op[tbox][s_idx], multipole[sbox][t_idx]);

          // Concat the leaf locals for the final result (L2T)
          // Permuted copy again
          int tb_idx = -1;
          for (target_body_type tb : bodies(tbox)) {
            ++tb_idx;
            result[tb.number()] += c[tb_idx];
          }
        }
      }
    }
  }

  // Check the result
  if (checkErrors) {
    std::cout << "Computing direct matvec..." << std::endl;

    // Compute the result with a direct matrix-vector multiplication
    std::vector<result_type> exact(M);
    Direct::matvec(K, sources, charges, targets, exact);

    double tot_error_sq = 0;
    double tot_norm_sq = 0;
    double tot_ind_rel_err = 0;
    double max_ind_rel_err = 0;
    for (unsigned k = 0; k < result.size(); ++k) {
      std::cout << result[k] << "\t" << exact[k] << std::endl;

      // Individual relative error
      double rel_error = norm_2(exact[k] - result[k]) / norm_2(exact[k]);
      tot_ind_rel_err += rel_error;
      // Maximum relative error
      max_ind_rel_err  = std::max(max_ind_rel_err, rel_error);

      // Total relative error
      tot_error_sq += norm_2_sq(exact[k] - result[k]);
      tot_norm_sq  += norm_2_sq(exact[k]);
    }
    double tot_rel_err = sqrt(tot_error_sq/tot_norm_sq);
    std::cout << "Vector  relative error: " << tot_rel_err << std::endl;

    double ave_rel_err = tot_ind_rel_err / result.size();
    std::cout << "Average relative error: " << ave_rel_err << std::endl;

    std::cout << "Maximum relative error: " << max_ind_rel_err << std::endl;
  }

  // TODO: Interpolative decompositions

  // TODO: Range based iterations for tree boxes -- simplification
  // TODO: Generalizations on Context so I can use it in this mock up.
}
