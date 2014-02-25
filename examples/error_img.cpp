#/** @file error.cpp
 * @brief Test the kernel and expansions by running an instance
 * of the kernel matrix-vector product and yielding statistics.
 */

#include <limits>       // std::numeric_limits

#include "fmmtl/KernelMatrix.hpp"
#include "fmmtl/Direct.hpp"

#include "UnitKernel.kern"
#include "ExpKernel.kern"

#include "LaplaceSpherical.hpp"
#include "LaplaceSpherical2.hpp"

#include "YukawaCartesian.hpp"

#define png_infopp_NULL (png_infopp)NULL
#define int_p_NULL (int*)NULL

#include <boost/gil/gil_all.hpp>
#include <boost/gil/extension/io/png_dynamic_io.hpp>

namespace gil = boost::gil;

/** Construct a rgb8_pixel_t from hue, saturation, and value.
 * @param[in] h hue (0 = red, 0.166 = yellow, 0.333 = green, 0.5 = cyan,
 *                   0.666 = blue, 0.833 = magenta)
 * @param[in] s saturation (0 = white, 1 = saturated)
 * @param[in] v value (0 = dark,  1 = bright)
 * @pre 0 <= @a h, @a s, @a v <= 1 */
gil::rgb8_pixel_t make_hsv(float h, float s, float v) {
  assert(0 <= h && h <= 1);
  assert(0 <= s && s <= 1);
  assert(0 <= v && v <= 1);

  if (s == 0)
    return gil::rgb8_pixel_t(255*v, 255*v, 255*v);

  float var_h = (h == 1 ? 0 : 6*h);
  int var_i = int(var_h);         // Or var_i = floor(var_h)
  float var_1 = v * (1 - s);
  float var_2 = v * (1 - s * (var_h - var_i));
  float var_3 = v * (1 - s * (1 - (var_h - var_i)));
  switch (var_i) {
    case 0: // hue [0,1/6): red <-> yellow
      return gil::rgb8_pixel_t(255*v, 255*var_3, 255*var_1);
    case 1: // hue [1/6,1/3): yellow <-> green
      return gil::rgb8_pixel_t(255*var_2, 255*v, 255*var_1);
    case 2: // hue [1/3,1/2): green <-> cyan
      return gil::rgb8_pixel_t(255*var_1, 255*v, 255*var_3);
    case 3: // hue [1/2,2/3): cyan <-> blue
      return gil::rgb8_pixel_t(255*var_1, 255*var_2, 255*v);
    case 4: // hue [2/3,5/6): blue <-> magenta
      return gil::rgb8_pixel_t(255*var_3, 255*var_1, 255*v);
    default: // hue [5/6,1): magenta <-> red
      return gil::rgb8_pixel_t(255*v, 255*var_1, 255*var_2);
  }
}

/** Construct a color suitable for heat map display.
 * @param[in] v Heat value (0 = cold, 1 = hot)
 * @pre 0 <= @a v <= 1
 *
 * The resulting color is a fully-saturated bright color ranging from
 * purple/blue for cold to red for hot. */
gil::rgb8_pixel_t make_heat(float v) {
  assert(0 <= v && v <= 1);
  return make_hsv(0.78f - 0.76f * v, 1, 1);
}



int main(int argc, char **argv)
{
  int N = 1 << 20;

  // Parse custom command line args
  for (int i = 1; i < argc; ++i) {
    if (strcmp(argv[i],"-N") == 0) {
      N = atoi(argv[++i]);
    }
  }
  // Round up to the nearest square number
  int n_side = int(std::ceil(std::sqrt(N)));
  N = n_side * n_side;

  // Init the FMM Kernel and options
  FMMOptions opts = get_options(argc, argv);
  //typedef UnitExpansion kernel_type;
  //typedef ExpExpansion kernel_type;
  typedef LaplaceSpherical2 kernel_type;
  //typedef YukawaCartesian kernel_type;

  // Init kernel
  kernel_type K;

  typedef kernel_type::point_type point_type;
  typedef kernel_type::source_type source_type;
  typedef kernel_type::target_type target_type;
  typedef kernel_type::charge_type charge_type;
  typedef kernel_type::result_type result_type;

  std::cout << "Initializing source and N = " << N << " targets..." << std::endl;

  // Init a square targets
  double xmin = -1;
  double xmax = 1;
  double ymin = -1;
  double ymax = 1;

  std::vector<target_type> targets;
  targets.reserve(N);
  for (int n = 0; n < n_side; ++n) {
    for (int m = 0; m < n_side; ++m) {
      targets.push_back(target_type(xmin + n * (xmax-xmin) / (n_side-1),
                                    ymin + m * (ymax-ymin) / (n_side-1)));
    }
  }
  //int middle = n_side/2 * n_side + n_side/2;
  int middle = fmmtl::random<unsigned>::get(0, N);

  // Init charges, only the middle source has a charge
  std::vector<charge_type> charges(N);
  charges[middle] = charge_type(1);

  std::cout << "Building the kernel matrix..." << std::endl;

  // Build the FMM
  fmmtl::kernel_matrix<kernel_type> A = K(targets, targets);
  A.set_options(opts);

  std::cout << "Performing the kernel matrix-vector mult..." << std::endl;

  // Execute the FMM
  std::vector<result_type> result = A * charges;

  // Check the result
  std::cout << "Computing direct kernel matrix-vector mult..." << std::endl;

  // Compute the result with a direct matrix-vector multiplication
  std::vector<result_type> exact(N);
  Direct::matvec(K,
                 targets.begin()+middle, targets.begin()+middle+1, charges.begin()+middle,
                 targets.begin(), targets.end(), exact.begin());

  std::cout << "Computing the errors..." << std::endl;

  std::vector<double> log_error(result.size());
  double min_error = std::numeric_limits<double>::max();
  double max_error = std::numeric_limits<double>::lowest();
  for (unsigned k = 0; k < result.size(); ++k) {
    double error = norm(exact[k] - result[k]) / norm(exact[k]);
    if (error > 1e-15)
      log_error[k] = std::log10(error);
    else
      log_error[k] = -16;
    min_error = std::min(min_error, log_error[k]);
    max_error = std::max(max_error, log_error[k]);
  }

  std::cout << "Min log error: " << min_error << std::endl;
  std::cout << "Max log error: " << max_error << std::endl;

  // Fill the image with the errors computed above
  gil::rgb8_image_t img(n_side, n_side);
  auto img_view = gil::view(img);
  for (int n = 0; n < n_side; ++n) {
    for (int m = 0; m < n_side; ++m) {
      img_view(n,m) = make_heat((log_error[n*n_side+m] - min_error) / (max_error - min_error));
    }
  }
  gil::png_write_view("fmmtl_errors.png", const_view(img));
}
