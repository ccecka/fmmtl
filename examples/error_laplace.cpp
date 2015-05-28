#/** @file error.cpp
 * @brief Test the kernel and expansions by running an instance
 * of the kernel matrix-vector product and yielding statistics.
 */

#include "fmmtl/KernelMatrix.hpp"
#include "fmmtl/Direct.hpp"
#include "fmmtl/util/Clock.hpp"
#include "fmmtl/numeric/random.hpp"

#include "LaplaceSpherical.hpp"
#include "LaplaceCartesian.hpp"


int main(int argc, char **argv)
{
  int N = 10000;
  bool checkErrors = true;

  // Parse custom command line args
  for (int i = 1; i < argc; ++i) {
    if (strcmp(argv[i],"-N") == 0) {
      N = atoi(argv[++i]);
    } else if (strcmp(argv[i],"-nocheck") == 0) {
      checkErrors = false;
    }
  }

  // Init the FMM Kernel and options
  FMMOptions opts = get_options(argc, argv);

  // Init kernel
  typedef LaplaceSpherical kernel_type;
  kernel_type K(5);
  //typedef LaplaceCartesian<5> kernel_type;
  //kernel_type K;

  typedef kernel_type::source_type source_type;
  typedef kernel_type::target_type target_type;
  typedef kernel_type::charge_type charge_type;
  typedef kernel_type::result_type result_type;

  // Init points and charges
  std::vector<source_type> points = fmmtl::random_n(N);

  std::vector<charge_type> charges = fmmtl::random_n(N);

  // Build the FMM
  fmmtl::kernel_matrix<kernel_type> A = K(points, points);
  A.set_options(opts);

  // Execute the FMM
  Clock t1;
  std::vector<result_type> result = A * charges;
  double time1 = t1.seconds();
  std::cout << "FMM in " << time1 << " secs" << std::endl;

  // Execute the FMM
  Clock t2;
  result = A * charges;
  double time2 = t2.seconds();
  std::cout << "FMM in " << time2 << " secs" << std::endl;

  // Execute the FMM
  Clock t3;
  result = A * charges;
  double time3 = t3.seconds();
  std::cout << "FMM in " << time3 << " secs" << std::endl;

  // Check the result
  if (checkErrors) {
    std::cout << "Computing direct matvec..." << std::endl;

    std::vector<result_type> exact(N);

    // Compute the result with a direct matrix-vector multiplication
    { ScopeClock timer("Direct in ");
    fmmtl::direct(K, points, charges, exact);
    }

    double tot_error_sq = 0;
    double tot_norm_sq = 0;
    double tot_ind_rel_err = 0;
    double max_ind_rel_err = 0;
    for (unsigned k = 0; k < result.size(); ++k) {
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
}
