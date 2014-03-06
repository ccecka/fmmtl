#include "fmmtl/KernelMatrix.hpp"
#include "fmmtl/Logger.hpp"

#include "LaplaceSpherical.hpp"

int main() {
  typedef LaplaceSpherical expansion_type;
  expansion_type K(5);
  typedef expansion_type::source_type source_type;
  typedef expansion_type::charge_type charge_type;
  typedef expansion_type::target_type target_type;
  typedef expansion_type::result_type result_type;

  FMMOptions opts;

  for (double n = 4; n <= 6; n += 0.125) {
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
    std::vector<result_type> r = A * charges;

    // Compute the product
    Ticker tick;
    std::vector<result_type> result = A * charges;
    double time = tick.seconds();

    std::cout << N << "\t" << time << "\n";
  }
}
