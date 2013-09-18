/** @file correctness.cpp
 * @brief Test the tree and tree traversal by running an instance
 * of the UnitKernel or ExpKernel with random points and charges
 */

#include "fmmtl/KernelMatrix.hpp"
#include "fmmtl/Direct.hpp"

#include "UnitKernel.kern"
#include "ExpKernel.kern"

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
  typedef ExpExpansion kernel_type;
  kernel_type K;

  typedef kernel_type::point_type point_type;
  typedef kernel_type::source_type source_type;
  typedef kernel_type::target_type target_type;
  typedef kernel_type::charge_type charge_type;
  typedef kernel_type::result_type result_type;

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

    int wrong_results = 0;
    for (unsigned k = 0; k < result.size(); ++k) {
      if ((exact[k] - result[k]) / exact[k] > 1e-10) {
        std::cout << "[" << std::setw(log10(N)+1) << k << "]"
                  << " Exact: " << exact[k]
                  << ", FMM: " << result[k] << std::endl;
        std::cout << (exact[k] - result[k]) / exact[k] << std::endl;
        ++wrong_results;
      }
    }
    std::cout << "Wrong counts: " << wrong_results << " of " << N << std::endl;
  }
}
