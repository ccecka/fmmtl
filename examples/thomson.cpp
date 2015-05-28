#/** @file thomson.cpp
 * @brief Use an FMM to optimize the total potential energy of N points
 * constrained to the surface of a 3-sphere.
 */

#include <iostream>
#include <fstream>

#include "fmmtl/KernelMatrix.hpp"
#include "fmmtl/Direct.hpp"
#include "fmmtl/util/Clock.hpp"
#include "fmmtl/numeric/random.hpp"

// FMM operator for the laplace potential/force (1/|x-y|)
#include "LaplaceSpherical.hpp"


int main(int argc, char **argv)
{
  // Number of points
  std::size_t N = 10000;

  // Parse custom command line args
  for (int i = 1; i < argc; ++i) {
    if (strcmp(argv[i],"-N") == 0) {
      N = atoi(argv[++i]);
    }
  }

  // Init the FMM Kernel and options
  FMMOptions opts = get_options(argc, argv);

  // Init kernel
  const unsigned p = 7;   // Order of the approximation
  typedef LaplaceSpherical kernel_type;
  kernel_type K(p);

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
      p = fmmtl::random<point_type>::get();  // point in [ 0,1]^DIM
      p = 2*p - 1;                           // point in [-1,1]^DIM
      mag = norm_2(p);
    } while (mag > 1);
    p /= mag;
  }
  // Init charges, could optimize out...
  std::vector<charge_type> charges(N, 1);


  std::cout << "Atoms Generated" << std::endl;

  long double lastE = 0;
  double scale = 1;

  for (unsigned iter = 1; iter < 300; ++iter) {
    // Build the FMM with the current points
    fmmtl::kernel_matrix<kernel_type> A = K(points, points);
    A.set_options(opts);

    // Compute forces and energy
    std::vector<result_type> pot_force = A * charges;


    // Compute the maximum cross product as a gradient magnitude proxy
    // and the total energy
    double max_cross_sq = 0;
    long double E = 0;
    for (std::size_t k = 0; k < N; ++k) {
      // Get the force part of the result
      const result_type& pf = pot_force[k];
      E += pf[0];

      point_type f = point_type(pf[1], pf[2], pf[3]);
      max_cross_sq = std::max(max_cross_sq, norm_2_sq(cross(points[k],f)));
    }
    // Full matvec double counts potential interactions
    E /= 2;

    // Step size heuristic for the gradient descent
    // TODO: Armijo/Wolfe/Raleigh
    if (lastE <= E)
      scale *= 0.5;
    else
      scale *= 1.1;
    lastE = E;

    // Use the scale and gradient magnitude proxy to compute the step size
    if (max_cross_sq == 0) break;
    double step = scale / (N * std::sqrt(max_cross_sq));

    // Print iteration data
    std::cout << "Iter: " << iter << "\n";
    std::cout << "  Grad:   " << std::setprecision(4) << max_cross_sq << "\n";
    std::cout << "  Scale:  " << std::setprecision(4) << scale << "\n";
    std::cout << "  Energy: " << std::setprecision(20) << E << std::endl;

    // Update from forces and normalize
    for (std::size_t k = 0; k < N; ++k) {
      for (std::size_t d = 0; d < 3; ++d)
        points[k][d] -= step * pot_force[k][1+d];
      points[k] /= norm_2(points[k]);
    }
  }

  // Write points to file
  std::ofstream file(std::string("thomson") + std::to_string(N));
  file << N << std::endl;
  for (auto&& p : points)
    file << std::setprecision(15) << p << std::endl;
}
