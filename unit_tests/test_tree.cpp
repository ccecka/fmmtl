#include "fmmtl/tree/Octree.hpp"
#include "fmmtl/Vec.hpp"

#include <iostream>
#include <string>

// Random number in [0,1)
inline double drand() {
  return ::drand48();
}


int main(int argc, char** argv)
{
  if (argc < 2) {
    std::cerr << "Usage: test_tree NUM_POINTS\n";
    exit(1);
  }

  int N = atoi(argv[1]);

  constexpr unsigned dimension = 2;
  typedef Vec<dimension,double> point_type;

  std::vector<point_type> points(N);
  for (int k = 0; k < N; ++k)
    for (unsigned d = 0; d < dimension; ++d)
      points[k][d] = drand();

  NDTree<dimension> tree(points.begin(), points.end());

  std::cout << tree << "\n";

  return 0;
}
