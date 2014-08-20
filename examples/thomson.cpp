#/** @file thomson.cpp
 * @brief Use an FMM to optimize the total potential energy of N points
 * constrained to the surface of a 3-sphere.
 */

#include <iostream>
#include <fstream>

#include "fmmtl/KernelMatrix.hpp"
#include "fmmtl/Direct.hpp"
#include "fmmtl/util/Clock.hpp"

// FMM operator for the laplace potential/force (1/|x-y|)
#include "LaplaceSpherical.hpp"


int main(int argc, char **argv)
{
  std::size_t N = 60000;
  bool checkErrors = true;

  // Parse custom command line args
  for (int i = 1; i < argc; ++i) {
    if (strcmp(argv[i],"-N") == 0) {
      N = atoi(argv[++i]);
    }
  }

  // Init the FMM Kernel and options
  FMMOptions opts = get_options(argc, argv);

  // Init kernel
  typedef LaplaceSpherical kernel_type;
  kernel_type K;

  typedef kernel_type::point_type  point_type;  // Vec<3,double>
  typedef kernel_type::charge_type charge_type; // double
  // Result type: (laplace potential, laplace force)
  typedef kernel_type::result_type result_type; // Vec<4,double>

  // Seed the RNG
  fmmtl::default_generator.seed(1337);

  // Init points uniformly on the sphere
  std::vector<point_type> points(N);
  for (auto&& p : points) {
    double mag;
    do {
      p = fmmtl::random<point_type>::get();
      mag = norm_2(p);
    } while (mag > 1);
    p /= mag;
  }
  // Init charges, could optimize out...
  std::vector<charge_type> charges(N, 1);


  std::cout << "Atoms Generated" << std::endl;

  long double E = 0;
  long double lastE = 0;
  double scale = 100;

  for (unsigned iter = 1; iter < 300; ++iter) {
    // Build the FMM with the current points
    fmmtl::kernel_matrix<kernel_type> A = K(points, points);
    A.set_options(opts);

    // Compute forces and energy
    std::vector<result_type> pot_force = A * charges;

    // Compute total energy
    E = 0;
    for (auto&& r : pot_force)
      E += r[0];

    // If the decrease is small, scale down
    // TODO: better scheduling/scaling for step size
    if (lastE - E <= 0)
      scale /= 2;
    lastE = E;
    scale *= 1.1;

    // Print iteration data
    std::cout << "Iter: " << iter << "\n";
    std::cout << "  Scale: " << scale << "\n";
    std::cout << "  Energy: " << std::setprecision(10) << E << std::endl;

    // Update from forces and normalize
    for (std::size_t k = 0; k < N; ++k) {
      for (std::size_t d = 0; d < points[k].size(); ++d)
        points[k][d] += scale * pot_force[k][1+d];
      points[k] /= norm_2(points[k]);
    }
  }

  // Write points to file
  std::ofstream file(std::string("thomson") + std::to_string(N));
  file << N << std::endl;
  for (auto&& p : points)
    file << std::setprecision(15) << p << std::endl;
}
