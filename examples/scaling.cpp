#include "fmmtl/KernelMatrix.hpp"

#include "LaplaceSpherical.hpp"

double get_time() {
  struct timeval tv;
  gettimeofday(&tv,nullptr);
  return (double)(tv.tv_sec+tv.tv_usec*1e-6);
}

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

  double tic, toc;

  for (double n = 4; n <= 6; n += 0.125) {
    int N = int(pow(10,n));

    std::vector<point_type> points(N);
    for (unsigned k = 0; k < points.size(); ++k)
      points[k] = point_type(drand(),drand(),drand());

    std::vector<charge_type> charges(N);
    for (unsigned k = 0; k < charges.size(); ++k)
      charges[k] = drand();

    // Initialize matrix
    fmm_matrix<expansion_type> plan = make_fmm_matrix(K, points, opts);

    // Compute the product
    tic = get_time();
    std::vector<result_type> result = plan * charges;
    toc = get_time();

    std::cout << N << "\t" << toc - tic << "\n";
  }
}
