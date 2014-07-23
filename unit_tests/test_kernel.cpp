#include "KernelSkeleton.kern"
#include "UnitKernel.kern"
#include "ExpKernel.kern"

#include "Laplace.kern"
#include "Yukawa.kern"
#include "Helmholtz.kern"
#include "Stokes.kern"
#include "Gaussian.kern"
#include "BiotSavart.kern"

#include "fmmtl/Direct.hpp"
#include "fmmtl/meta/kernel_traits.hpp"

template <class Kernel>
void test_kernel(const Kernel& K) {
  typedef typename Kernel::source_type source_type;
  typedef typename Kernel::charge_type charge_type;
  typedef typename Kernel::target_type target_type;
  typedef typename Kernel::result_type result_type;

  result_type r = K(target_type(), source_type()) * charge_type();

  std::cout << std::endl;
  std::cout << typeid(Kernel).name() << ":" << std::endl;
  std::cout << r << std::endl;
  std::cout << KernelTraits<Kernel>() << std::endl;
}



int main() {
  test_kernel(KernelSkeleton()); // KernelSkeleton
  test_kernel(UnitPotential());  // UnitKernel
  test_kernel(ExpPotential());   // ExpKernel

  test_kernel(LaplacePotential()); // Laplace
  test_kernel(LaplaceKernel());    // Laplace

  test_kernel(YukawaPotential()); // Yukawa
  test_kernel(YukawaKernel());    // Yukawa

  test_kernel(HelmholtzPotential()); // Helmholtz
  test_kernel(HelmholtzKernel());    // Helmholtz

  test_kernel(Stokeslet()); // Stokes

  test_kernel(Gaussian<1>()); // Gaussian
  test_kernel(Gaussian<4>()); // Gaussian

  test_kernel(BiotSavart());      // BiotSavart
  test_kernel(RosenheadMoore());  // RosenheadMoore
}
