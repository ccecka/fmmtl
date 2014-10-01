#include <iostream>
#include <string>

#include "fmmtl/tree/KDTree3.hpp"
#include "fmmtl/numeric/Vec.hpp"

#include "fmmtl/util/Clock.hpp"

int main(int argc, char** argv)
{
  if (argc < 2) {
    std::cerr << "Usage: " << argv[0] << " NUM_POINTS\n";
    exit(1);
  }

  int N = atoi(argv[1]);

  const std::size_t dim = 3;

  typedef Vec<dim,double> point_type;

  std::vector<point_type> points(N);
  for (int k = 0; k < N; ++k)
    points[k] = fmmtl::random<point_type>::get();

  Clock timer;

  timer.start();
  fmmtl::KDTree<dim> tree(points, 100);
  double time = timer.seconds();
  std::cout << "KDTree Construction: " << time << std::endl;

  std::cout << tree.bounding_box() << std::endl;
  std::cout << tree << "\n";
}
