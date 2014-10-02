#include <iostream>
#include <string>

#include "fmmtl/numeric/Vec.hpp"
#include "fmmtl/tree/KDTree3.hpp"

#include "fmmtl/util/Clock.hpp"

int main(int argc, char** argv)
{
 if (argc < 2) {
    std::cerr << "Usage: " << argv[0] << " [-N NUM_POINTS] [-no-print]\n";
    exit(1);
  }

  int N = 10000;
  bool print = true;

  // Parse custom command line args
  for (int i = 1; i < argc; ++i) {
    if (strcmp(argv[i],"-N") == 0) {
      N = atoi(argv[++i]);
    } else if (strcmp(argv[i],"-no-print") == 0) {
      print = false;
    }
  }

  const std::size_t dim = 3;

  typedef Vec<dim,double> point_type;

  std::vector<point_type> points = fmmtl::random_n(N);

  Clock timer;
  fmmtl::KDTree<dim> tree(points, 10);
  double time = timer.seconds();
  std::cout << "KDTree Construction: " << time << std::endl;

  if (print) {
    std::cout << tree.bounding_box() << std::endl;
    std::cout << tree << "\n";
  }
}
