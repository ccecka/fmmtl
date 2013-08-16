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

  template <class Context, class BoxPairIter>
  P2P_Compressed<typename Context::kernel_type>
  make(Context& c, BoxPairIter first, BoxPairIter last) {
    /*
    typename Context::source_iterator first_source = bc.source_begin();
    typename Context::source_iterator first_target = bc.target_begin();

    unsigned num_targets = bc.target_tree().bodies();
    //unsigned num_sources = bc.source_tree().bodies();

    // Interaction list for each target box
    // target_first -> {(source_first, source_last), ...}
    // TODO: faster?
    typedef std::pair<unsigned, unsigned> upair;
    std::vector<std::vector<upair>> target2sources(num_targets);
    // A list of target ranges we've seen: {(target_first, target_last), ...}
    std::vector<upair> target_ranges;

    for ( ; first != last; ++first) {
      const box_pair& b2b = *first;
      const source_box_type& source_box = b2b.first;
      const target_box_type& target_box = b2b.second;

			// Target range
			unsigned i_begin = bc.target_begin(target_box) - first_target;

      auto& sources = target2sources[i_begin];
      if (sources.empty()) {
        // First time we've seen this target range
        unsigned i_end = bc.target_end(target_box) - first_target;
        target_ranges.push_back(upair(i_begin, i_end));
      }
      //FMMTL_ASSERT(targets.find(upair(i_begin, bc.target_end(target_box)-first_target)) != targets.end());

			// Source range
			unsigned j_begin = bc.source_begin(source_box) - first_source;
			unsigned j_end = bc.source_end(source_box) - first_source;
      sources.push_back(upair(j_begin,j_end));
    }

    // Construct a compressed interaction list
    std::vector<unsigned> target_ptr(target_ranges.size() + 1);
    auto target_ptr_curr = target_ptr.begin();

    std::vector<upair> source_ranges(p2p_list.size());
    auto source_ranges_curr = source_ranges.begin();

    // For all the target ranges
    for (auto& target_range : target_ranges) {
      // Record the offset for this source range
      *target_ptr_curr = source_ranges_curr - source_ranges.begin();
      ++target_ptr_curr;

      // Copy the interacting source ranges
      auto& sources = target2sources[target_range.first];
      source_ranges_curr = std::copy(sources.begin(), sources.end(),
                                     source_ranges_curr);
    }

    *target_ptr_curr = source_ranges_curr - source_ranges.begin();

    // Sanity checking
    FMMTL_ASSERT(*target_ptr_curr == source_ranges.size());
    FMMTL_ASSERT(++target_ptr_curr == target_ptr.end());

    // TODO
		return make_p2p_gpu(bc.kernel(),
                        target_ranges,
                        target_ptr,
                        source_ranges,
												first_source, bc.source_end(),
												first_target, bc.target_end());
    */
    return P2P_Compressed<typename Context::kernel_type>();
  }
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
