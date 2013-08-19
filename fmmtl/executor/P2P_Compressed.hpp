#pragma once
/** @file EvalP2P_GPU.hpp
 * @brief Header file for the GPU_P2P class.
 *
 * Note: This header file is compiled with nvcc and must use C++03.
 */

#include <vector>

#include <iostream>

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
  void* data_;

  // Device data for P2P computation
  std::pair<unsigned,unsigned>* target_ranges_;
  unsigned* source_range_ptrs_;
  std::pair<unsigned,unsigned>* source_ranges_;

  // Device source and target arrays
  source_type* sources_;
  target_type* targets_;

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
  static
  P2P_Compressed<typename Context::kernel_type>*
  make(Context& c, BoxPairIter first, BoxPairIter last) {
    typename Context::source_iterator first_source = c.source_begin();
    typename Context::source_iterator first_target = c.target_begin();

    unsigned num_targets = c.target_tree().bodies();
    //unsigned num_sources = c.source_tree().bodies();
    unsigned num_box_pairs = last - first;

    // Interaction list for each target box
    // (target_first,target_last) -> {(source_first, source_last), ...}
    // TODO: faster?
    typedef std::pair<unsigned, unsigned> upair;
    std::vector<std::vector<upair> > target2sources(num_targets);
    // A list of target ranges we've seen: {(target_first, target_last), ...}
    std::vector<upair> target_ranges;

    for ( ; first != last; ++first) {
      const typename BoxPairIter::value_type& bpair = *first;
      const typename Context::source_box_type& source_box = bpair.first;
      const typename Context::target_box_type& target_box = bpair.second;

      // Target boxes need to be leaf boxes because the computations are
      // grouped by disjoint target ranges
      // TODO: Generalize?
      FMMTL_ASSERT(target_box.is_leaf());

			// Target range
			unsigned i_begin = c.target_begin(target_box) - first_target;
      unsigned i_end   = c.target_end(target_box) - first_target;

			// Source range
			unsigned j_begin = c.source_begin(source_box) - first_source;
			unsigned j_end   = c.source_end(source_box) - first_source;

      // If this is the first time we've seen this target range
      if (target2sources[i_begin].empty())
        target_ranges.push_back(upair(i_begin, i_end));

      target2sources[i_begin].push_back(upair(j_begin,j_end));
    }

    unsigned num_target_ranges = target_ranges.size();

    // Construct a compressed interaction list
    std::vector<unsigned> target_ptr(num_target_ranges + 1);
    target_ptr[0] = 0;
    std::vector<upair> source_ranges(num_box_pairs);
    std::vector<upair>::iterator source_ranges_curr = source_ranges.begin();

    // For all the target ranges
    for (unsigned k = 0; k < num_target_ranges; ++k) {
      // Copy the source ranges that interact with the kth target range
      unsigned i_begin = target_ranges[k].first;
      source_ranges_curr = std::copy(target2sources[i_begin].begin(),
                                     target2sources[i_begin].end(),
                                     source_ranges_curr);

      // Record the stop index
      target_ptr[k+1] = source_ranges_curr - source_ranges.begin();
    }

    // Sanity checking
    FMMTL_ASSERT(target_ptr.back() == source_ranges.size());
    FMMTL_ASSERT(source_ranges_curr == source_ranges.end());

    // Copy the source and target ranges into contiguous vectors
    std::vector<source_type> sources(c.source_begin(), c.source_end());
    std::vector<target_type> targets(c.target_begin(), c.target_end());

    return new P2P_Compressed<Kernel>(target_ranges,
                                      target_ptr,
                                      source_ranges,
                                      sources,
                                      targets);
  }
};
