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

#include "fmmtl/Direct.hpp"

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
  int N = 5000;
  int M = 5000;
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
  std::vector<result_type> result(M, result_type());

  // Dimension of the tree
  const unsigned D = fmmtl::dimension<source_type>::value;
  static_assert(D == fmmtl::dimension<target_type>::value, "Dimension mismatch");

  // Construct two trees
  fmmtl::NDTree<D> source_tree(sources, 2);
  fmmtl::NDTree<D> target_tree(targets, 2);

  typedef typename fmmtl::NDTree<D>::box_type target_box_type;
  typedef typename fmmtl::NDTree<D>::box_type source_box_type;

  // Operators associated with the source and target boxes
  std::vector<std::vector<matrix<double>>> p_op(source_tree.boxes());
  std::vector<matrix<kernel_value_type>> b_op(target_tree.boxes());

  //
  // Butterfly precomputation
  //

  // For all the levels of the target_tree
  //int max_level = std::min(source_tree.levels(), target_tree.levels()) - 1;
  int max_level = 3;

  std::cout << "max_level = " << max_level << std::endl;
  for (int level = 0; level <= max_level; ++level) {

    std::cout << "Level " << level << std::endl;

    // For all boxes in the opposing level of the source tree
    int slevel = max_level - level;
    auto sb_end = source_tree.box_end(slevel);
    for (auto sbi = source_tree.box_begin(slevel); sbi != sb_end; ++sbi) {
      source_box_type sbox = *sbi;

      p_op[sbox.index()].resize(std::distance(target_tree.box_begin(level),
                                              target_tree.box_end(level)));

      assert(sbox.num_children() == (1 << D));

      // For all the boxes in this level of the target tree
      auto tb_end = target_tree.box_end(level);
      unsigned tb_offset = (*(target_tree.box_begin(level))).index();
      for (auto tbi = target_tree.box_begin(level); tbi != tb_end; ++tbi) {
        target_box_type tbox = *tbi;

        assert(tbox.num_children() == (1 << D));

        std::cout << "Precomputing " << level << "  " << tbox.index() << "  " << sbox.index() << std::endl;

        // Interpolative decomposition of tbox x sbox
        // or propogation of previous interpolative decompositions

        if (level == 0) {
          // Form interpolative decomposition (dummy)
          p_op[sbox.index()][tbox.index()-tb_offset] = identity_matrix<double>(sbox.num_bodies(), sbox.num_bodies());
          //t_op[tbox.index()] = B-matrix... ?
        } else {
          // Get the "B"-matrix of the target's parent and partition according to child partitions
          // Form the interpolative decomposition of each partition (dummy)
          p_op[sbox.index()][tbox.index()-tb_offset] = identity_matrix<double>(sbox.num_bodies(), sbox.num_bodies());

          // Since we've only done dummy P-matrices, this will be the full evaluated rows
          if (level == max_level) {
            // Construct the elements of the kernel matrix
            // Could use a ublas kernel_matrix class here...
            auto& k = b_op[tbox.index()] = matrix<kernel_value_type>(tbox.num_bodies(), sbox.num_bodies());
            std::cout << "Matrix " << k.size1() << "," << k.size2() << std::endl;
            unsigned t_offset = (*(tbox.body_begin())).index();
            unsigned s_offset = (*(sbox.body_begin())).index();
            for (auto tb = tbox.body_begin(); tb != tbox.body_end(); ++tb) {
              target_type& t = targets[(*tb).number()];
              for (auto sb = sbox.body_begin(); sb != sbox.body_end(); ++sb) {
                source_type& s = sources[(*sb).number()];
                k((*tb).index() - t_offset, (*sb).index() - s_offset) = K(t, s);
              }
            }
          }
        }
      }
    }
  }

  //
  // Butterfly Application
  //

  std::vector<std::vector<vector<charge_type>>> multipole(source_tree.boxes());
  std::vector<vector<kernel_value_type>> local(target_tree.boxes());

  // For all the levels of the target_tree
  for (int level = 0; level <= max_level; ++level) {

    // For all boxes in the opposing level of the source tree
    int slevel = max_level - level;
    auto sb_end = source_tree.box_end(slevel);
    for (auto sbi = source_tree.box_begin(slevel); sbi != sb_end; ++sbi) {
      source_box_type sbox = *sbi;

      multipole[sbox.index()].resize(p_op[sbox.index()].size());

      // For all the boxes in this level of the target tree
      auto tb_end = target_tree.box_end(level);
      unsigned tb_offset = (*(target_tree.box_begin(level))).index();
      for (auto tbi = target_tree.box_begin(level); tbi != tb_end; ++tbi) {
        target_box_type tbox = *tbi;

        std::cout << "Computing " << level << "  " << tbox.index() << "  " << sbox.index() << std::endl;

        // Application of precomputed interpolative decompositions
        if (level == 0) {
          // Create the "multipole" from the sources (S2M)
          // Construct the local source vector explictly...
          vector<charge_type> c(sbox.num_bodies());
          unsigned s_offset = (*(sbox.body_begin())).index();
          for (auto si = sbox.body_begin(); si != sbox.body_end(); ++si) {
            auto sbody = *si;
            c[sbody.index()-s_offset] = charges[sbody.number()];
          }
          // Compute the product (dummy)
          multipole[sbox.index()][tbox.index()-tb_offset] = prod(p_op[sbox.index()][tbox.index()-tb_offset], c);
        } else {
          // Concatenate the "multipoles" from the children and apply the op (M2M)
          vector<charge_type> c(sbox.num_bodies());
          unsigned s_offset = (*(sbox.body_begin())).index();
          unsigned tbp_offset = (*(target_tree.box_begin(level-1))).index();
          for (auto ci = sbox.child_begin(); ci != sbox.child_end(); ++ci) {
            auto cbox = *ci;
            unsigned c_offset = (*(cbox.body_begin())).index();
            for (auto cb = cbox.body_begin(); cb != cbox.body_end(); ++cb) {
              auto cbody = *cb;
              c[cbody.index()-s_offset] = multipole[cbox.index()][tbox.parent().index()-tbp_offset][cbody.index()-c_offset];
            }
          }
          // And take the product
          multipole[sbox.index()][tbox.index()-tb_offset] = prod(p_op[sbox.index()][tbox.index()-tb_offset], c);
        }

        if (level == max_level) {
          // Take the product of the ID and the "multipole" (M2L)
          local[tbox.index()] = prod(b_op[tbox.index()], multipole[sbox.index()][tbox.index()-tb_offset]);
        }
      }
    }
  }

  // Concatenate the leaf locals for the final result (L2T)

  // For all the boxes in this level of the target tree
  auto tb_end = target_tree.box_end(max_level);
  for (auto tbi = target_tree.box_begin(max_level); tbi != tb_end; ++tbi) {
    target_box_type tbox = *tbi;

    // For all the bodies
    unsigned t_offset = (*(tbox.body_begin())).index();
    for (auto tb = tbox.body_begin(); tb != tbox.body_end(); ++tb) {
      auto b = *tb;
      result[b.number()] = local[tbox.index()][b.index() - t_offset];
    }
  }

  // Check the result
  if (checkErrors) {
    std::cout << "Computing direct matvec..." << std::endl;

    // Compute the result with a direct matrix-vector multiplication
    std::vector<result_type> exact(M, result_type());
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
