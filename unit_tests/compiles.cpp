#include "fmmtl/KernelMatrix.hpp"

#include "KernelSkeleton.kern"

// Random number in [0,1)
inline double drand() {
  return ::drand48();
}

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

  typedef KernelExpansion expansion_type;
  typedef expansion_type::point_type point_type;
  typedef expansion_type::source_type source_type;
  typedef expansion_type::target_type target_type;
  typedef expansion_type::charge_type charge_type;
  typedef expansion_type::result_type result_type;

  expansion_type K;
  FMMOptions opts = get_options(argc, argv);

  std::vector<source_type> sources(N);
  for (source_type& s : sources)
    s = source_type(drand(), drand(), drand());

  std::vector<target_type> targets(M);
  for (target_type& t : targets)
    t = target_type(drand(), drand(), drand());

  std::vector<charge_type> charges(N);
  std::vector<result_type> result;

  // Single Tree
  fmm_matrix<expansion_type> A1 = K(sources, sources);
  result = A1 * charges;
  //for (result_type r : result) std::cout << r << std::endl;

  // Dual Tree
  fmm_matrix<expansion_type> A2  = K(targets, sources);
  result = A2 * charges;
  //for (result_type r : result) std::cout << r << std::endl;
}
