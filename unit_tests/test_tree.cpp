#include "tree/Octree.hpp"
#include "Vec.hpp"

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

  typedef Vec<3,double> point_type;

  std::vector<point_type> points(N);
  for (int k = 0; k < N; ++k)
    points[k] = point_type(drand(), drand(), drand());

  Octree<point_type> otree(points.begin(), points.end());

  std::cout << otree << "\n";

  return 0;
}
