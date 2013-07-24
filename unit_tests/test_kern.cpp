#include "KernelSkeleton.kern"
#include "UnitKernel.kern"
#include "ExpKernel.kern"

#include "Laplace.kern"
#include "Yukawa.kern"
#include "Helmholtz.kern"
#include "Stokes.kern"

#include "fmmtl/Direct.hpp"
#include "fmmtl/KernelTraits.hpp"

template <class Kernel>
void test_kernel(const Kernel& K) {
  typedef typename Kernel::source_type source_type;
  typedef typename Kernel::charge_type charge_type;
  typedef typename Kernel::target_type target_type;
  typedef typename Kernel::result_type result_type;

  std::vector<source_type> sources(1);
  std::vector<charge_type> charges(1);
  std::vector<target_type> targets(1);
  std::vector<result_type> results(1);

  Direct::matvec(K, sources, charges, targets, results);

  std::cout << results[0] << std::endl;
  std::cout << KernelTraits<Kernel>() << std::endl;
  std::cout << typeid(Kernel).name() << " OK." << std::endl;
}



int main() {
  test_kernel(Kernel());        // KernelSkeleton
  test_kernel(UnitPotential()); // UnitKernel
  test_kernel(ExpPotential());  // ExpKernel

  test_kernel(LaplacePotential()); // Laplace
  test_kernel(LaplaceKernel());    // Laplace

  test_kernel(YukawaPotential(1.0)); // Yukawa
  test_kernel(YukawaKernel(1.0));    // Yukawa

  test_kernel(HelmholtzPotential(1.0)); // Helmholtz
  test_kernel(HelmholtzKernel(1.0));    // Helmholtz

  test_kernel(Stokeslet()); // Stokes
}
