#include "fmmtl/KernelMatrix.hpp"
#include "LaplaceSpherical.hpp"

#include "fmmtl/Logger.hpp"

// Random number in [0,1)
inline double drand() {
  return ::drand48();
}

int main() {
  typedef LaplaceSpherical expansion_type;
  expansion_type K(5);
  typedef expansion_type::point_type point_type;
  typedef expansion_type::charge_type charge_type;
  typedef expansion_type::result_type result_type;

  FMMOptions opts;

  for (double n = 4; n <= 6; n += 0.125) {
    int N = int(pow(10,n));

    std::vector<point_type> points(N);
    for (unsigned k = 0; k < points.size(); ++k)
      points[k] = point_type(drand(),drand(),drand());

    std::vector<charge_type> charges(N);
    for (unsigned k = 0; k < charges.size(); ++k)
      charges[k] = drand();

    // Initialize matrix
    fmm_matrix<expansion_type> A = make_fmm_matrix(K, points, opts);

    // Compute the product
    Ticker tick;
    std::vector<result_type> result = A * charges;
    double time = tick.seconds();

    std::cout << N << "\t" << time << "\n";
  }
}
