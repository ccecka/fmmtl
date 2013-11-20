#include "fmmtl/tree/NDTree.hpp"
#include "fmmtl/Vec.hpp"

#include <iostream>
#include <string>

int main(int argc, char** argv)
{
  if (argc < 2) {
    std::cerr << "Usage: " << argv[0] << " NUM_POINTS\n";
    exit(1);
  }

  int N = atoi(argv[1]);

  typedef Vec<2,double> point_type;

  std::vector<point_type> points(N);
  for (int k = 0; k < N; ++k)
    points[k] = fmmtl::random<point_type>::get();

  NDTree<2> tree(points.begin(), points.end());

  std::cout << tree << "\n";
}
