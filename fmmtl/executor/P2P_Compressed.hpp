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

  // Supporting data
  void* data;

  // Device data for P2P computation
  unsigned* target_ranges;
  unsigned* source_range_ptrs;
  unsigned* source_ranges;

  // Device source and target arrays
  source_type* sources;
  target_type* targets;

  P2P_Compressed();

	P2P_Compressed(std::vector<std::pair<unsigned,unsigned> >& target_ranges,
                 std::vector<unsigned>& target_ptrs,
                 std::vector<std::pair<unsigned,unsigned> >& source_ranges,
                 std::vector<source_type>& sources,
                 std::vector<target_type>& targets);

  ~P2P_Compressed();

  template <class Context>
  void execute(Context& c) {
    return execute(c.kernel(),
                   c.charge_begin(), c.charge_end(),
                   c.result_begin(), c.result_end());
  }

  template <typename ChargeIter, typename ResultIter>
        void execute(const Kernel& K,
                     ChargeIter cfirst, ChargeIter clast,
                     ResultIter rfirst, ResultIter rlast) {
    // TODO: Ugh, iterator type hiding via copy
    std::vector<charge_type> charges(cfirst, clast);
    std::vector<result_type> results(rfirst, rlast);

    execute(K, charges, results);

    std::copy(results.begin(), results.end(), rfirst);
  }

  void execute(const Kernel& K,
               const std::vector<charge_type>& charges,
               std::vector<result_type>& results);

  static void execute(const Kernel& K,
                      const std::vector<source_type>& s,
                      const std::vector<charge_type>& c,
                      const std::vector<target_type>& t,
                      std::vector<result_type>& r);
};



template <typename Kernel, typename SourceIter, typename TargetIter>
P2P_Compressed<Kernel>
make_p2p_gpu(const Kernel&,
             std::vector<std::pair<unsigned,unsigned> >& target_ranges,
             std::vector<unsigned>& target_ptrs,
             std::vector<std::pair<unsigned,unsigned> >& source_ranges,
             SourceIter sfirst, SourceIter slast,
             TargetIter tfirst, TargetIter tlast) {
  // TODO: Ugh, iterator type hiding via copy
  std::vector<typename Kernel::source_type> sources(sfirst, slast);
  std::vector<typename Kernel::target_type> targets(tfirst, tlast);

	return P2P_Compressed<Kernel>(target_ranges,
                                target_ptrs,
                                source_ranges,
                                sources,
                                targets);
}
