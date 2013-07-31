/** @file correctness.cpp
 * @brief Test the tree and tree traversal by running an instance
 * of the UnitKernel with random points and charges
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
  int N = 10000;  // num sources
  int M = 10000;  // num targets
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

  // Init sources
  std::vector<source_type> sources(N);
  for (unsigned k = 0; k < sources.size(); ++k)
    sources[k] = source_type(drand(), drand(), drand());
  // Init charges
  std::vector<charge_type> charges(N);
  for (unsigned k = 0; k < charges.size(); ++k)
    charges[k] = drand();

  // Init targets
  std::vector<target_type> targets(M);
  for (unsigned k = 0; k < targets.size(); ++k)
    targets[k] = target_type(drand(), drand(), drand());

  // Build the FMM
  fmm_matrix<kernel_type> A = make_fmm_matrix(K, sources, targets, opts);

  // Execute the FMM
  std::vector<result_type> result = A * charges;

  // Check the result
  if (checkErrors) {
    std::cout << "Computing direct matvec..." << std::endl;

    std::vector<result_type> exact(M);

    // Compute the result with a direct matrix-vector multiplication
    Direct::matvec(K, sources, charges, targets, exact);

    int wrong_results = 0;
    for (unsigned k = 0; k < result.size(); ++k) {
      if ((exact[k] - result[k]) / exact[k] > 1e-10) {
        std::cout << "[" << std::setw(log10(M)+1) << k << "]"
                  << " Exact: " << exact[k]
                  << ", FMM: " << result[k] << std::endl;
        std::cout << (exact[k] - result[k]) / exact[k] << std::endl;
        ++wrong_results;
      }
    }
    std::cout << "Wrong counts: " << wrong_results << " of " << M << std::endl;
  }
}
