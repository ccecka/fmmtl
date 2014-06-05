#include "fmmtl/numeric/Vec.hpp"

#include <iostream>

int main() {
  std::cout << "Is POD: " << std::is_pod<Vec<3,double> >::value << std::endl;
  std::cout << "Is trivial: " << std::is_trivial<Vec<3,double> >::value << std::endl;
  std::cout << "Is standard layout: " << std::is_standard_layout<Vec<3,double> >::value << std::endl;

  Vec<3,double> v0 = {};
  std::cout << v0 << std::endl;

  Vec<3,double> v1 = {0., 1., 2};

  Vec<3,double> v2 = {3, 4, 5};

  Vec<3,double> v3 = v1 * (v2 + v2);

  std::cout << v3 << std::endl;
  std::cout << norm_2(v3) << std::endl;

  Vec<4, double> v = {1, 2.1f, 3.14, 2u};

  std::cout << norm_1(v) << std::endl;
  std::cout << norm_2(v) << std::endl;
  std::cout << norm_inf(v) << std::endl;
  std::cout << v << std::endl;
}
