#include "fmmtl/tree/NDTree.hpp"
#include "fmmtl/numeric/Vec.hpp"

#include <iostream>
#include <string>

int main(int argc, char** argv)
{
  if (argc < 2) {
    std::cerr << "Usage: " << argv[0] << " NUM_POINTS\n";
    exit(1);
  }

  int N = atoi(argv[1]);

  typedef Vec<3,double> point_type;

  std::vector<point_type> points = fmmtl::random_n(N);

  fmmtl::NDTree<3> tree(points.begin(), points.end());

  std::cout << tree.bounding_box() << std::endl;
  std::cout << tree << "\n";
}
