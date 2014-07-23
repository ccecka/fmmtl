#include <iostream>

#include <thrust/device_vector.h>

#include "fmmtl/numeric/Vec.hpp"

int main() {
  typedef Vec<3,double> source_type;

  thrust::device_vector<source_type> d_s(10);

  for (unsigned k = 0; k < d_s.size(); ++k)
    std::cout << (source_type) d_s[k] << std::endl;
}
