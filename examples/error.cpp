#/** @file error.cpp
 * @brief Test the kernel and expansions by running an instance
 * of the kernel matrix-vector product and yielding statistics.
 */

#include "fmmtl/KernelMatrix.hpp"
#include "fmmtl/Direct.hpp"

#include "UnitKernel.kern"
#include "ExpKernel.kern"

#include "LaplaceSpherical.hpp"
#include "YukawaCartesian.hpp"

#include "Math.hpp"

// Random number in [0,1)
inline double drand() {
  return ::drand48();
}

// Random number in [A,B)
inline double drand(double A, double B) {
  return (B-A) * drand() + A;
}


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
  //typedef UnitExpansion kernel_type;
  //typedef ExpExpansion kernel_type;
  //typedef LaplaceSpherical kernel_type;
  typedef YukawaCartesian kernel_type;

  typedef kernel_type::point_type point_type;
  typedef kernel_type::source_type source_type;
  typedef kernel_type::target_type target_type;
  typedef kernel_type::charge_type charge_type;
  typedef kernel_type::result_type result_type;

  // Init kernel
  kernel_type K;

  // Init points and charges
  std::vector<source_type> points(N);
  for (unsigned k = 0; k < points.size(); ++k)
    points[k] = source_type(drand(-1,1), drand(-1,1), drand(-1,1));

  std::vector<charge_type> charges(N);
  for (unsigned k = 0; k < charges.size(); ++k)
    charges[k] = drand(-1,1);

  // Build the FMM
  fmmtl::kernel_matrix<kernel_type> A = K(points, points);
  A.set_options(opts);

  // Execute the FMM
  std::vector<result_type> result = A * charges;

  // Check the result
  if (checkErrors) {
    std::cout << "Computing direct matvec..." << std::endl;

    std::vector<result_type> exact(N);

    // Compute the result with a direct matrix-vector multiplication
    Direct::matvec(K, points, charges, exact);

    double tot_error_sq = 0;
    double tot_norm_sq = 0;
    double tot_ind_rel_err = 0;
    double max_ind_rel_err = 0;
    for (unsigned k = 0; k < result.size(); ++k) {
      // Individual relative error
      double rel_error = norm(exact[k] - result[k]) / norm(exact[k]);
      tot_ind_rel_err += rel_error;
      // Maximum relative error
      max_ind_rel_err  = std::max(max_ind_rel_err, rel_error);

      // Total relative error
      tot_error_sq += normSq(exact[k] - result[k]);
      tot_norm_sq  += normSq(exact[k]);
    }
    double tot_rel_err = sqrt(tot_error_sq/tot_norm_sq);
    std::cout << "Vector  relative error: " << tot_rel_err << std::endl;

    double ave_rel_err = tot_ind_rel_err / result.size();
    std::cout << "Average relative error: " << ave_rel_err << std::endl;

    std::cout << "Maximum relative error: " << max_ind_rel_err << std::endl;
  }
}
