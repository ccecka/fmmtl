#include "fmmtl/KernelMatrix.hpp"

#include "KernelSkeleton.kern"

int main(int argc, char** argv) {
  int N = 100;
  int M = 100;
  // Parse custom command line args
  for (int i = 1; i < argc; ++i) {
    if (strcmp(argv[i],"-N") == 0) {
      N = atoi(argv[++i]);
    } else if (strcmp(argv[i],"-M") == 0) {
      M = atoi(argv[++i]);
    }
  }

  typedef ExpansionSkeleton expansion_type;
  typedef expansion_type::source_type source_type;
  typedef expansion_type::target_type target_type;
  typedef expansion_type::charge_type charge_type;
  typedef expansion_type::result_type result_type;

  expansion_type K;
  FMMOptions opts = get_options(argc, argv);

  std::vector<source_type> sources(N);
  for (source_type& s : sources)
    s = fmmtl::random<source_type>::get();

  std::vector<target_type> targets(M);
  for (target_type& t : targets)
    t = fmmtl::random<target_type>::get();

  std::vector<charge_type> charges(N);
  std::vector<result_type> result;

  // Single Tree
  fmmtl::kernel_matrix<expansion_type> A1 = K(sources, sources);
  A1.set_options(opts);
  result = A1 * charges;
  //for (result_type r : result) std::cout << r << std::endl;

  // Dual Tree
  fmmtl::kernel_matrix<expansion_type> A2  = K(targets, sources);
  A2.set_options(opts);
  result = A2 * charges;
  //for (result_type r : result) std::cout << r << std::endl;
}
