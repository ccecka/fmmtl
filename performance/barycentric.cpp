#include <iostream>
#include <iomanip>

#include "fmmtl/KernelMatrix.hpp"
#include "fmmtl/Direct.hpp"
#include "fmmtl/util/Clock.hpp"
#include "fmmtl/numeric/random.hpp"

#include "BarycentricTaylor.hpp"

int main(int argc, char** argv) {
  // The BarycentricTaylor expansion with polynomial order 20
  typedef BarycentricTaylor<20> expansion_type;
  expansion_type K;
  typedef expansion_type::source_type source_type;
  typedef expansion_type::charge_type charge_type;
  typedef expansion_type::target_type target_type;
  typedef expansion_type::result_type result_type;

  FMMOptions opts = get_options(argc, argv);

  auto tab = std::setw(15);
  std::cout << tab
            << "N" << tab
            << "Theta" << tab
            << "NCrit" << tab
            << "MaxRelError" << tab
            << "First" << tab
            << "Average" << tab
            << "Direct" << std::endl;

  for (double n = 1; n <= 5; n += 0.125) {
    int N = int(pow(10,n));

    std::vector<source_type> points = fmmtl::random_n(N);

    std::vector<charge_type> charges = fmmtl::random_n(N);

    // Initialize matrix
    fmmtl::kernel_matrix<expansion_type> A = K(points, points);
    A.set_options(opts);

    // Warm up
    Clock clock;

    clock.start();
    std::vector<result_type> result = A * charges;
    double init_time = clock.seconds();

    // Compute the product
    double time = 0;
    const unsigned ITER = 10;
    for (unsigned iter = 0; iter < ITER; ++iter) {
      clock.start();
      std::vector<result_type> result = A * charges;
      time += clock.seconds();
    }
    time /= ITER;

    // Compute the result with a direct matrix-vector multiplication
    clock.start();
    std::vector<result_type> exact(N);
    fmmtl::direct(K, points, charges, exact);
    double direct_time = clock.seconds();

    // Compute the error
    double max_ind_rel_err = 0;
    for (unsigned k = 0; k < result.size(); ++k) {
      // Individual relative error
      double rel_error = norm_2(exact[k] - result[k]) / norm_2(exact[k]);
      // Maximum relative error
      max_ind_rel_err  = std::max(max_ind_rel_err, rel_error);
    }

    std::cout << tab
              << N << tab
              << opts.theta << tab
              << opts.ncrit << tab
              << max_ind_rel_err << tab
              << init_time << tab
              << time << tab
              << direct_time << std::endl;
  }

  return 0;
}
