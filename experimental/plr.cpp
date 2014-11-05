
#include "fmmtl/numeric/Vec.hpp"
#include "fmmtl/numeric/norm.hpp"

#include "plr.hpp"

#include "fmmtl/Direct.hpp"

#include "fmmtl/util/Clock.hpp"


struct MyKernel {
  typedef double value_type;

  typedef Vec<2,value_type> source_type;
  typedef Vec<2,value_type> target_type;
  typedef value_type        charge_type;

  //typedef double result_type;
  //typedef double kernel_value_type;
  typedef std::complex<double>   result_type;
  typedef std::complex<double>   kernel_value_type;

  kernel_value_type operator()(const target_type& t,
                               const source_type& s) const {
    return std::exp(-0.5*norm_2_sq(s - t));
  }
};


int main(int argc, char** argv) {
  unsigned N = 1000;
  unsigned M = 1000;
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

  // PLR Parameters for the submatrix blocks
  static constexpr std::size_t max_rank = 5;
  static constexpr double      eps_tol  = 1e-3;

  // Define a kernel to work with
  using kernel_type = MyKernel;
  kernel_type kernel;

  using source_type = typename kernel_type::source_type;
  using target_type = typename kernel_type::target_type;
  using charge_type = typename kernel_type::charge_type;
  using result_type = typename kernel_type::result_type;

  using kernel_value_type = typename kernel_type::kernel_value_type;

  // Construct a set of sources and targets
  std::vector<source_type> sources = fmmtl::random_n(M);
  std::vector<target_type> targets = fmmtl::random_n(N);

  // Construct a set of charges and results
  std::vector<charge_type> charges = fmmtl::random_n(M);
  std::vector<result_type> results(N);

  // Get the dimension of the sources and targets
  constexpr unsigned SD = fmmtl::dimension<source_type>::value;
  constexpr unsigned TD = fmmtl::dimension<target_type>::value;

  // Generate the matrix
  std::vector<kernel_value_type> A(N*M);
  for (unsigned i = 0; i < N; ++i)
    for (unsigned j = 0; j < M; ++j)
      A[i*M+j] = kernel(targets[i], sources[j]);

  // Call the PLR Compression
  auto plr_plan
      = plr_compression<TD,SD>(A.data(), N, M,
                               reinterpret_cast<double*>(sources.data()),
                               reinterpret_cast<double*>(targets.data()),
                               max_rank, eps_tol);

  // Perform the matvec
  { ScopeClock timer("PLR MatVec: ");

    prod_acc(plr_plan, charges, results);

  } // timer


  // Check the result
  if (checkErrors) {
    std::cout << "Computing direct matvec..." << std::endl;

    // Compute the result with a direct matrix-vector multiplication
    std::vector<result_type> exact(M);
    { ScopeClock timer("Direct: ");
      fmmtl::direct(kernel, sources, charges, targets, exact);
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
    double tot_rel_err = std::sqrt(tot_error_sq/tot_norm_sq);
    std::cout << "Vector  relative error: " << tot_rel_err << std::endl;

    double ave_rel_err = tot_ind_rel_err / results.size();
    std::cout << "Average relative error: " << ave_rel_err << std::endl;

    std::cout << "Maximum relative error: " << max_ind_rel_err << std::endl;
  }
}
