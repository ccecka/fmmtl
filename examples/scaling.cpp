#include "fmmtl/KernelMatrix.hpp"
#include "fmmtl/util/Clock.hpp"

#include "LaplaceSpherical.hpp"

int main() {
  typedef LaplaceSpherical expansion_type;
  expansion_type K(5);
  typedef expansion_type::source_type source_type;
  typedef expansion_type::charge_type charge_type;
  typedef expansion_type::target_type target_type;
  typedef expansion_type::result_type result_type;

  FMMOptions opts;

  for (double n = 4; n <= 7; n += 0.125) {
    int N = int(pow(10,n));

    std::vector<source_type> points(N);
    for (unsigned k = 0; k < points.size(); ++k)
      points[k] = fmmtl::random<source_type>::get();

    std::vector<charge_type> charges(N);
    for (unsigned k = 0; k < charges.size(); ++k)
      charges[k] = fmmtl::random<charge_type>::get();

    // Initialize matrix
    fmmtl::kernel_matrix<expansion_type> A = K(points, points);
    A.set_options(opts);

    // Warm up
    double time = 0;
    Clock clock;

    clock.start();
    std::vector<result_type> r = A * charges;
    double init_time = clock.seconds();

    // Compute the product
    const unsigned ITER = 10;
    for (unsigned iter = 0; iter < ITER; ++iter) {
      clock.start();
      std::vector<result_type> result = A * charges;
      time += clock.seconds();
    }
    time /= ITER;

    std::cout << N << "\t" << init_time << "\t" << time << "\n";
  }

  return 0;
}
