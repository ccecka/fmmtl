#pragma once
/** @file EvalP2P_GPU.hpp
 * @brief Header file for the GPU_P2P class.
 *
 * Note: This header file is compiled with nvcc and must use C++03.
 */

#include <vector>

template <typename Kernel>
class P2P_Compressed {
 public:
  typedef Kernel kernel_type;

  // No KernelTraits here... only C++03
  typedef typename kernel_type::source_type source_type;
  typedef typename kernel_type::target_type target_type;
  typedef typename kernel_type::charge_type charge_type;
  typedef typename kernel_type::result_type result_type;

  // The kernel to apply in the P2P
  const kernel_type& K_;

  // Device data for P2P computation
  unsigned* d_target_ranges;
  unsigned* d_target_ptrs;
  unsigned* d_source_ranges;

  // Device source and target arrays
  source_type* d_sources;
  target_type* d_targets;

  // Supporting data
  void* data;

  P2P_Compressed(const Kernel& K);

	P2P_Compressed(const Kernel& K,
                 std::vector<std::pair<unsigned,unsigned> >& target_ranges,
                 std::vector<unsigned>& target_ptrs,
                 std::vector<std::pair<unsigned,unsigned> >& source_ranges,
                 std::vector<source_type>& sources,
                 std::vector<target_type>& targets);

  template <typename ChargeIter, typename ResultIter>
  void execute(ChargeIter cfirst, ChargeIter clast,
               ResultIter rfirst, ResultIter rlast) {
    // TODO: Ugh, iterator type hiding via copy
    std::vector<charge_type> charges(cfirst, clast);
    std::vector<result_type> results(rfirst, rlast);

    execute(&charges[0], &results[0]);

    std::copy(results.begin(), results.end(), rfirst);
  }

  void execute(const charge_type* charges,
               result_type* results);
};



template <typename Expansion, typename SourceIter, typename TargetIter>
P2P_Compressed<typename Expansion::kernel_type>
make_p2p_gpu(const Expansion& K,
             std::vector<std::pair<unsigned,unsigned> >& target_ranges,
             std::vector<unsigned>& target_ptrs,
             std::vector<std::pair<unsigned,unsigned> >& source_ranges,
             SourceIter sfirst, SourceIter slast,
             TargetIter tfirst, TargetIter tlast) {
  // TODO: Ugh, iterator type hiding via copy
  std::vector<typename Expansion::source_type> sources(sfirst, slast);
  std::vector<typename Expansion::target_type> targets(tfirst, tlast);

	return P2P_Compressed<typename Expansion::kernel_type>(K,
                                                         target_ranges,
                                                         target_ptrs,
                                                         source_ranges,
                                                         sources,
                                                         targets);
}
