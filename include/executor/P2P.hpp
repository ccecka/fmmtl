#pragma once
/** @file P2P.hpp
 * @brief Dispatch methods for the P2P stage
 *
 */

#include "KernelTraits.hpp"
#include <type_traits>

struct P2P
{
  /** If no other P2P dispatcher matches */
  template <typename Kernel, typename... Args>
  inline static void eval(const Kernel&, Args...) {
    std::cerr << "Kernel does not have a correct op() or P2P!\n";
    std::cerr << KernelTraits<Kernel>() << std::endl;
    exit(1);
  }

  /** Dual-Evaluation dispatch
   */
  template <typename Kernel>
  inline static
  typename std::enable_if<KernelTraits<Kernel>::has_eval_op &
                          !KernelTraits<Kernel>::has_transpose>::type
  eval(const Kernel& K,
       const typename Kernel::source_type& p1,
       const typename Kernel::charge_type& c1,
       typename Kernel::result_type& r1,
       const typename Kernel::source_type& p2,
       const typename Kernel::charge_type& c2,
       typename Kernel::result_type& r2)
  {
    r1 += K(p1,p2) * c2;
    r2 += K(p2,p1) * c1;
  }

  /** Dual-Evaluation dispatch
   */
  template <typename Kernel>
  inline static
  typename std::enable_if<KernelTraits<Kernel>::has_eval_op &
                          KernelTraits<Kernel>::has_transpose>::type
  eval(const Kernel& K,
       const typename Kernel::source_type& p1,
       const typename Kernel::charge_type& c1,
       typename Kernel::result_type& r1,
       const typename Kernel::source_type& p2,
       const typename Kernel::charge_type& c2,
       typename Kernel::result_type& r2)
  {
    typedef typename Kernel::kernel_value_type kernel_value_type;

    kernel_value_type k12 = K(p1,p2);
    r1 += k12 * c2;
    kernel_value_type k21 = K.transpose(k12);
    r2 += k21 * c1;
  }

	/** Asymmetric vectorized P2P dispatch
	 */
	template <typename Kernel,
	          typename SourceIter, typename ChargeIter,
	          typename TargetIter, typename ResultIter>
	inline static
	typename std::enable_if<KernelTraits<Kernel>::has_vector_P2P_asymm>::type
	eval(const Kernel& K,
       SourceIter s_first, SourceIter s_last, ChargeIter c_first,
       TargetIter t_first, TargetIter t_last, ResultIter r_first)
	{
		K.P2P(s_first, s_last, c_first,
		      t_first, t_last, r_first);
	}

  /** Symmetric vectorized P2P dispatch
   * @pre source_type == target_type
   */
  template <typename Kernel,
            typename SourceIter, typename ChargeIter, typename ResultIter>
  inline static
  typename std::enable_if<KernelTraits<Kernel>::has_vector_P2P_symm>::type
  eval(const Kernel& K,
       SourceIter p1_first, SourceIter p1_last, ChargeIter c1_first,
       ResultIter r1_first,
       SourceIter p2_first, SourceIter p2_last, ChargeIter c2_first,
       ResultIter r2_first)
  {
    K.P2P(p1_first, p1_last, c1_first,
          p2_first, p2_last, c2_first,
          r1_first, r2_first);
  }

	/** Asymmetric P2P using the evaluation operator
	 * r_i += sum_j K(t_i, s_j) * c_j
	 *
	 * @param[in] ...
	 */
	template <typename Kernel,
	          typename SourceIter, typename ChargeIter,
	          typename TargetIter, typename ResultIter>
	inline static
  typename std::enable_if<KernelTraits<Kernel>::has_eval_op &
                          !KernelTraits<Kernel>::has_vector_P2P_asymm>::type
  eval(const Kernel& K,
       SourceIter s_first, SourceIter s_last, ChargeIter c_first,
       TargetIter t_first, TargetIter t_last, ResultIter r_first)
  {
    typedef typename Kernel::source_type source_type;
    typedef typename Kernel::charge_type charge_type;
    typedef typename Kernel::target_type target_type;
    typedef typename Kernel::result_type result_type;

    // TODO?
    // Optimize on if(std::iterator_traits<All Iters>::iterator_category == random_access_iterator)
    // to eliminate multiple increments

    static_assert(std::is_same<target_type,
                               typename TargetIter::value_type>::value,
                  "TargetIter::value_type != Kernel::target_type");
    static_assert(std::is_same<source_type,
                               typename SourceIter::value_type>::value,
                  "SourceIter::value_type != Kernel::source_type");
    static_assert(std::is_same<charge_type,
                               typename ChargeIter::value_type>::value,
                  "ChargeIter::value_type != Kernel::charge_type");
    static_assert(std::is_same<result_type,
                               typename ResultIter::value_type>::value,
                  "ResultIter::value_type != Kernel::result_type");

    for ( ; t_first != t_last; ++t_first, ++r_first) {
      const target_type& t = *t_first;
      result_type& r       = *r_first;

      SourceIter s = s_first;
      ChargeIter c = c_first;
      for ( ; s != s_last; ++s, ++c)
        r += K(t,*s) * (*c);
    }
  }

  /** Symmetric P2P, off-diagonal block using the evaluation operator
   * r2_i += sum_j K(p2_i, p1_j) * c1_j
   * r1_j += sum_i K(p1_j, p2_i) * c2_i
   *
   * @param[in] ...
   * @pre source_type == target_type
   * @pre For all i,j we have p1_i != p2_j
   */
  template <typename Kernel,
            typename SourceIter, typename ChargeIter, typename ResultIter>
  inline static
  typename std::enable_if<!KernelTraits<Kernel>::has_vector_P2P_symm>::type
  eval(const Kernel& K,
       SourceIter p1_first, SourceIter p1_last, ChargeIter c1_first,
       ResultIter r1_first,
       SourceIter p2_first, SourceIter p2_last, ChargeIter c2_first,
       ResultIter r2_first)
  {
    typedef typename Kernel::source_type source_type;
    typedef typename Kernel::charge_type charge_type;
    typedef typename Kernel::target_type target_type;
    typedef typename Kernel::result_type result_type;

    static_assert(std::is_same<source_type,target_type>::value,
                  "source_type != target_type in symmetric P2P");
    static_assert(std::is_same<source_type,
                               typename SourceIter::value_type>::value,
                  "SourceIter::value_type != Kernel::source_type");
    static_assert(std::is_same<charge_type,
                               typename ChargeIter::value_type>::value,
                  "ChargeIter::value_type != Kernel::charge_type");
    static_assert(std::is_same<result_type,
                               typename ResultIter::value_type>::value,
                  "ResultIter::value_type != Kernel::result_type");

    // TODO
    // Optimize on random_access_iterator?

    for ( ; p1_first != p1_last; ++p1_first, ++c1_first, ++r1_first) {
      const source_type& p1 = *p1_first;
      const charge_type& c1 = *c1_first;
      result_type& r1       = *r1_first;

      SourceIter p2i = p2_first;
      ChargeIter c2i = c2_first;
      ResultIter r2i = r2_first;
      for ( ; p2i != p2_last; ++p2i, ++c2i, ++r2i)
        P2P::eval(K, p1, c1, r1, *p2i, *c2i, *r2i);
    }
  }

  /** Symmetric P2P, diagonal block using the evaluation operator
   * r_i += sum_j K(p_i, p_j) * c_j
   *
   * @pre source_type == target_type
   */
  template <typename Kernel,
            typename SourceIter, typename ChargeIter, typename ResultIter>
  inline static
  typename std::enable_if<!KernelTraits<Kernel>::has_vector_P2P_symm>::type
  eval(const Kernel& K,
       SourceIter p_first, SourceIter p_last,
       ChargeIter c_first, ResultIter r_first)
  {
    typedef typename Kernel::source_type source_type;
    typedef typename Kernel::charge_type charge_type;
    typedef typename Kernel::target_type target_type;
    typedef typename Kernel::result_type result_type;

    static_assert(std::is_same<source_type, target_type>::value,
                  "source_type != target_type in symmetric P2P");
    static_assert(std::is_same<source_type,
                               typename SourceIter::value_type>::value,
                  "SourceIter::value_type != Kernel::source_type");
    static_assert(std::is_same<charge_type,
                               typename ChargeIter::value_type>::value,
                  "ChargeIter::value_type != Kernel::charge_type");
    static_assert(std::is_same<result_type,
                               typename ResultIter::value_type>::value,
                  "ResultIter::value_type != Kernel::result_type");

    // TODO
    // Optimize on random_access_iterator?

    SourceIter pi = p_first;
    ChargeIter ci = c_first;
    ResultIter ri = r_first;
    for ( ; pi != p_last; ++pi, ++ci, ++ri) {
      const source_type& p = *pi;
      const charge_type& c = *ci;
      result_type& r       = *ri;

      // The off-diagonal elements
      SourceIter pj = p_first;
      ChargeIter cj = c_first;
      ResultIter rj = r_first;
      for ( ; pj != pi; ++pj, ++cj, ++rj)
        P2P::eval(K, p, c, r, *pj, *cj, *rj);

      // The diagonal element
      r += K(p,p) * c;
    }
  }


  //////////////////////////////////////
  /////// Context Dispatchers //////////
  //////////////////////////////////////

  struct ONE_SIDED {};
  struct TWO_SIDED {};

  /** Asymmetric P2P
   */
  template <typename Context>
  inline static void eval(Context& c,
                          const typename Context::source_box_type& source,
                          const typename Context::target_box_type& target,
                          const ONE_SIDED&)
  {
#ifdef DEBUG
    std::cout << "P2P:\n  " << source << "\n  " << target << std::endl;
#endif

    P2P::eval(c.kernel(),
              c.source_begin(source), c.source_end(source),
              c.charge_begin(source),
              c.target_begin(target), c.target_end(target),
              c.result_begin(target));
  }

  /** Symmetric P2P
   */
  template <typename Context>
  inline static void eval(Context& c,
                          const typename Context::source_box_type& box1,
                          const typename Context::target_box_type& box2,
                          const TWO_SIDED&)
  {
#ifdef DEBUG
    std::cout << "P2P:\n  " << box1 << "\n  " << box2 << std::endl;
    std::cout << "P2P:\n  " << box2 << "\n  " << box1 << std::endl;
#endif

    P2P::eval(c.kernel(),
              c.source_begin(box1), c.source_end(box1),
              c.charge_begin(box1), c.result_begin(box1),
              c.target_begin(box2), c.target_end(box2),
              c.charge_begin(box2), c.result_begin(box2));
  }

  /** Symmetric P2P
   */
  template <typename Context>
  inline static void eval(Context& c,
                          const typename Context::source_box_type& box)
  {
#ifdef DEBUG
    std::cout << "P2P:\n  " << box << std::endl;
#endif

    P2P::eval(c.kernel(),
              c.source_begin(box), c.source_end(box),
              c.charge_begin(box), c.result_begin(box),
              c.target_begin(box), c.target_end(box),
              c.charge_begin(box), c.result_begin(box));
  }
};



/**
 * Batched P2P methods
 **/

#include <cmath>

#include "Evaluator.hpp"
#include "P2P_Compressed.hpp"

#define BOOST_UBLAS_NDEBUG
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/vector.hpp>
namespace ublas = boost::numeric::ublas;

/** A lazy P2P evaluator which saves a list of pairs of boxes
 * That are sent to the P2P dispatcher on demand.
 */
template <typename Context>
class P2P_Batch
	: public EvaluatorBase<Context> {
  //! Kernel type
  typedef typename Context::kernel_type kernel_type;
  //! Kernel value type
  typedef typename Context::kernel_value_type kernel_value_type;

  //! Type of box
  typedef typename Context::box_type box_type;
  //! Box list for P2P interactions    TODO: could further compress these...
  typedef std::pair<box_type,box_type> box_pair;
  std::vector<box_pair> p2p_list;

 public:
  /** Insert a source-target box interaction to the interaction list */
  void insert(const box_type& box1, const box_type& box2) {
    p2p_list.push_back(std::make_pair(box1,box2));
  }

  /** Compute all interations in the interaction list */
  virtual void execute(Context& bc) const {
    for (const box_pair& b2b : p2p_list)
      P2P::eval(bc.kernel(), bc, b2b.first, b2b.second, P2P::ONE_SIDED());
  }

  class P2P_Matrix
      : public EvaluatorBase<Context> {
    ublas::compressed_matrix<kernel_value_type> A;

   public:
    virtual void execute(Context& bc) {
      // printf("EvalLocalSparse::execute(Context&)\n");

      typedef typename Context::charge_type charge_type;
      ublas::vector<charge_type> charges(bc.source_tree().bodies());
      std::copy(bc.charge_begin(), bc.charge_end(), charges.begin());

      // Call the matvec
      typedef typename Context::result_type result_type;
      ublas::vector<result_type> results = ublas::prod(A, charges);

      // Accumulate results
      std::transform(results.begin(), results.end(),
                     bc.result_begin(), bc.result_begin(),
                     std::plus<result_type>());
    }
  };


  /** Convert the interaction list to an interaction matrix
   * by evaluating all of the elements.
   */
  ublas::compressed_matrix<kernel_value_type> to_matrix(Context& bc) {
    auto first_source = bc.source_begin();
    auto first_target = bc.target_begin();

    // Interaction list for each target body
    unsigned rows = bc.target_tree().bodies();
    unsigned cols = 0;
    unsigned nnz = 0;
    std::vector<std::vector<unsigned>> csr(rows);

    for (const box_pair& b2b : p2p_list) {
      const box_type& box1 = b2b.first;
      const box_type& box2 = b2b.second;

      auto source_end = bc.source_end(box1);
      auto target_end = bc.target_end(box2);
      for (auto t = bc.target_begin(box2); t != target_end; ++t) {
        // Row
        unsigned i = t - first_target;
        std::vector<unsigned>& csri = csr[i];

        for (auto s = bc.source_begin(box1); s != source_end; ++s) {
          // Column
          unsigned j = s - first_source;

          //assert(std::find(csri.begin(), csri.end(), j) == csri.end());
          csri.push_back(j);
          ++nnz;
          cols = std::max(cols, j);
        }
      }
    }
    ++cols;

    // The precomputed interaction matrix
    ublas::compressed_matrix<kernel_value_type> m(rows, cols, nnz);

    typedef typename kernel_type::source_type source_type;
    typedef typename kernel_type::target_type target_type;
    for (unsigned i = 0; i < csr.size(); ++i) {
      // Insert elements to compressed_matrix in-order for efficiency
      std::sort(csr[i].begin(), csr[i].end());
      const target_type& target = first_target[i];

      for (unsigned j : csr[i]) {
        const source_type& source = first_source[j];
        m.push_back(i, j, bc.kernel()(target, source));
      }
    }

    return m;
  }

	/** All boxes interactions have been inserted, stage for GPU P2P
	 */
  P2P_Compressed<typename kernel_type::kernel_type> to_gpu(Context& bc) {
    auto first_source = bc.source_begin();
    auto first_target = bc.target_begin();

    unsigned num_targets = bc.target_tree().bodies();
    //unsigned num_sources = bc.source_tree().bodies();

    // Interaction list for each target box
    // target_first -> {(source_first, source_last), ...}
    // TODO: faster?
    typedef std::pair<unsigned, unsigned> upair;
    std::vector<std::vector<upair>> target2sources(num_targets);
    // A list of target ranges we've seen: {(target_first, target_last), ...}
    std::vector<upair> target_ranges;

    for (const box_pair& b2b : p2p_list) {
      const box_type& source_box = b2b.first;
      const box_type& target_box = b2b.second;

			// Target range
			unsigned i_begin = bc.target_begin(target_box) - first_target;

      auto& sources = target2sources[i_begin];
      if (sources.empty()) {
        // First time we've seen this target range
        unsigned i_end = bc.target_end(target_box) - first_target;
        target_ranges.push_back(upair(i_begin, i_end));
      }
      //assert(targets.find(upair(i_begin, bc.target_end(target_box)-first_target)) != targets.end());

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
    assert(*target_ptr_curr == source_ranges.size());
    assert(++target_ptr_curr == target_ptr.end());

    // TODO
		return make_p2p_gpu(bc.kernel(),
                        target_ranges,
                        target_ptr,
                        source_ranges,
												first_source, bc.source_end(),
												first_target, bc.target_end());
  }
};
