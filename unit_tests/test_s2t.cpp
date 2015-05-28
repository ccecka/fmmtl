#include "KernelSkeleton.kern"
#include "UnitKernel.kern"
#include "ExpKernel.kern"

#include "Laplace.kern"
#include "Yukawa.kern"
#include "Helmholtz.kern"
#include "Stokes.kern"

#include "fmmtl/dispatch/S2T.hpp"

#include "fmmtl/numeric/random.hpp"

#include <vector>
#include <set>

typedef std::pair<unsigned,unsigned> upair;


/** Create random ranges inside the range [min,max) such that
 * the resulting ranges are of maximum size ncrit
 *
 * @param min,max  The minimum and maximum of the original range: [min,max)
 * @param ncrit    The maximum size of any resulting range:
 *                 for all i, result[i].second - result[i].first < ncrit
 */
std::vector<upair> random_ranges(unsigned min, unsigned max, unsigned ncrit) {
  std::set<unsigned> idx;
  idx.insert(min);
  idx.insert(max);
  auto it = idx.begin();
  while (*it != max) {
    unsigned curr = *it;
    auto cit = it;
    unsigned next = *(++cit);
    if (next - curr >= ncrit)
      idx.insert(curr + fmmtl::random<unsigned>::get() % (next-curr));
    else
      ++it;
  }

  std::vector<upair> ranges;
  auto i = idx.begin();
  while (*i != max) {
    upair r;
    r.first  = *i;
    r.second = *(++i);
    ranges.push_back(r);
  }

  return ranges;
}

template <typename T1, typename T2>
void print_results(const std::vector<T1>& exact, const std::vector<T2>& result) {
  double tot_error_sq = 0;
  double tot_norm_sq = 0;
  double tot_ind_rel_err = 0;
  double max_ind_rel_err = 0;
  for (unsigned k = 0; k < result.size(); ++k) {
    // Individual relative error
    double rel_error = norm_2(exact[k] - result[k]) / norm_2(exact[k]);
    //if (rel_error > 1e-14) std::cout << k << std::endl;

    tot_ind_rel_err += rel_error;
    // Maximum relative error
    max_ind_rel_err  = std::max(max_ind_rel_err, rel_error);

    // Total relative error
    tot_error_sq += norm_2_sq(exact[k] - result[k]);
    tot_norm_sq  += norm_2_sq(exact[k]);
  }
  double tot_rel_err = sqrt(tot_error_sq/tot_norm_sq);
  std::cout << "  Vector  relative error: " << tot_rel_err << std::endl;

  double ave_rel_err = tot_ind_rel_err / result.size();
  std::cout << "  Average relative error: " << ave_rel_err << std::endl;

  std::cout << "  Maximum relative error: " << max_ind_rel_err << std::endl;
}



int main(int argc, char** argv) {
  unsigned N = 10000;           // N rows (# targets)
  unsigned M = 10000;           // M cols (# sources)
  unsigned range_size = 256;    // Maximum block size

  // Parse custom command line args
  for (int i = 1; i < argc; ++i) {
    if (strcmp(argv[i],"-N") == 0) {
      N = atoi(argv[++i]);
    } else if (strcmp(argv[i],"-M") == 0) {
      M = atoi(argv[++i]);
    }
  }

  // Define the Kernel function we'll be testing
  typedef LaplaceKernel kernel_type;
  // Create a Kernel
  kernel_type K;

  // Get the types from the Kernel
  typedef typename kernel_type::source_type source_type;
  typedef typename kernel_type::charge_type charge_type;
  typedef typename kernel_type::target_type target_type;
  typedef typename kernel_type::result_type result_type;

  // Create a vector of sources
  std::vector<source_type> s = fmmtl::random_n(M);

  // Create a vector of charges
  std::vector<charge_type> c = fmmtl::random_n(M);

  // Create a vector of targets
  std::vector<target_type> t = fmmtl::random_n(N);

  //**************************************//

  // Create a vector of results
  std::vector<result_type> r_cpu(N);

  // Compute the full S2T on the CPU.
  fmmtl::detail::block_eval(K,
                            s.begin(), s.end(), c.begin(),
                            t.begin(), t.end(), r_cpu.begin());

#if defined(FMMTL_DEBUG)
  std::cout << "CPU:" << std::endl;
  for (result_type& ri : r_cpu) std::cout << ri << std::endl;
#endif

  //**************************************//

  // Create a vector of results
  std::vector<result_type> r_gpu(N);

  // Compute the full S2T on the GPU
  S2T_Compressed<kernel_type>::execute(K, s, c, t, r_gpu);

#if defined(FMMTL_DEBUG)
  std::cout << "GPU:" << std::endl;
  for (result_type& ri : r_gpu) std::cout << ri << std::endl;
#endif

  //**************************************//

  // Chop [0,N) and [0,M) into ranges
  std::vector<upair> source_range = random_ranges(0, s.size(), range_size);
  std::vector<upair> target_range = random_ranges(0, t.size(), range_size);

#if defined(FMMTL_DEBUG)
  for (upair& pi : source_range)
    std::cout << pi.first << ", " << pi.second << std::endl;
  for (upair& pi : target_range)
    std::cout << pi.first << ", " << pi.second << std::endl;
#endif

  std::vector<upair> sr_list(source_range.size() * target_range.size());
  std::vector<upair> tr_list(source_range.size() * target_range.size());

  // Create all pairs of source_ranges and target_ranges for testing
  for (unsigned si = 0; si != source_range.size(); ++si) {
    for (unsigned ti = 0; ti != target_range.size(); ++ti) {
      sr_list[si * target_range.size() + ti] = source_range[si];
      tr_list[si * target_range.size() + ti] = target_range[ti];
    }
  }

#if defined(FMMTL_DEBUG)
  for (unsigned k = 0; k < sr_list.size(); ++k) {
    std::cout << sr_list[k].first << ", " << sr_list[k].second << " :: "
              << tr_list[k].first << ", " << tr_list[k].second << std::endl;
  }
#endif

  // Create the blocked S2T on the GPU
  S2T_Compressed<kernel_type>* p2p
      = S2T_Compressed<kernel_type>::make(sr_list.begin(), sr_list.end(),
                                          tr_list.begin(),
                                          s, t);
  // Create a vector of results
  std::vector<result_type> r_gpu_b(N);

  // Execute the blocked S2T on the GPU
  p2p->execute(K, c, r_gpu_b);

#if defined(FMMTL_DEBUG)
  std::cout << "GPU Blocked:" << std::endl;
  for (result_type& ri : r_gpu_b) std::cout << ri << std::endl;
#endif

  // Print results
  std::cout << "CPU-GPU:" << std::endl;
  print_results(r_cpu, r_gpu);
  std::cout << std::endl;

  std::cout << "CPU-GPU Blocked:" << std::endl;
  print_results(r_cpu, r_gpu_b);
}
