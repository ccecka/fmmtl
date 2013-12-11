#pragma once
/** @file P2P.hpp
 * @brief Dispatch methods for the P2P stage
 *
 */

#include "fmmtl/Logger.hpp"
#include "fmmtl/meta/kernel_traits.hpp"
#include <iterator>
#include <type_traits>

struct P2P
{
  /** Dual-Evaluation dispatch
   */
  template <typename Kernel>
  inline static
  typename std::enable_if<KernelTraits<Kernel>::has_eval_op &
                          !KernelTraits<Kernel>::has_transpose>::type
  symm_eval(const Kernel& K,
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
  symm_eval(const Kernel& K,
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
    r2 += K.transpose(k12) * c1;
  }

	/** Asymmetric block P2P dispatch
	 */
	template <typename Kernel,
	          typename SourceIter, typename ChargeIter,
	          typename TargetIter, typename ResultIter>
	inline static
	typename std::enable_if<KernelTraits<Kernel>::has_vector_P2P_asymm>::type
	block_eval(const Kernel& K,
             SourceIter s_first, SourceIter s_last, ChargeIter c_first,
             TargetIter t_first, TargetIter t_last, ResultIter r_first)
	{
		K.P2P(s_first, s_last, c_first,
		      t_first, t_last, r_first);
	}

	/** Asymmetric block P2P using the evaluation operator
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
  block_eval(const Kernel& K,
             SourceIter s_first, SourceIter s_last, ChargeIter c_first,
             TargetIter t_first, TargetIter t_last, ResultIter r_first)
  {
    typedef typename Kernel::source_type source_type;
    typedef typename Kernel::charge_type charge_type;
    typedef typename Kernel::target_type target_type;
    typedef typename Kernel::result_type result_type;

    static_assert(std::is_same<source_type,
                  typename std::iterator_traits<SourceIter>::value_type
                  >::value, "SourceIter::value_type != Kernel::source_type");
    static_assert(std::is_same<charge_type,
                  typename std::iterator_traits<ChargeIter>::value_type
                  >::value, "ChargeIter::value_type != Kernel::charge_type");
    static_assert(std::is_same<target_type,
                  typename std::iterator_traits<TargetIter>::value_type
                  >::value, "TargetIter::value_type != Kernel::target_type");
    static_assert(std::is_same<result_type,
                  typename std::iterator_traits<ResultIter>::value_type
                  >::value, "ResultIter::value_type != Kernel::result_type");

    for ( ; t_first != t_last; ++t_first, ++r_first) {
      const target_type& t = *t_first;
      result_type& r       = *r_first;

      SourceIter si = s_first;
      ChargeIter ci = c_first;
      for ( ; si != s_last; ++si, ++ci)
        r += K(t,*si) * (*ci);
    }
  }

  /** Symmetric off-diagonal block P2P dispatch
   * @pre source_type == target_type
   */
  template <typename Kernel,
            typename SourceIter, typename ChargeIter, typename ResultIter>
  inline static
  typename std::enable_if<KernelTraits<Kernel>::has_vector_P2P_symm>::type
  block_eval(const Kernel& K,
             SourceIter p1_first, SourceIter p1_last, ChargeIter c1_first,
             ResultIter r1_first,
             SourceIter p2_first, SourceIter p2_last, ChargeIter c2_first,
             ResultIter r2_first)
  {
    K.P2P(p1_first, p1_last, c1_first,
          p2_first, p2_last, c2_first,
          r1_first, r2_first);
  }

  /** Symmetric off-diagonal block P2P using the evaluation operator
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
  block_eval(const Kernel& K,
             SourceIter p1_first, SourceIter p1_last, ChargeIter c1_first,
             ResultIter r1_first,
             SourceIter p2_first, SourceIter p2_last, ChargeIter c2_first,
             ResultIter r2_first)
  {
    typedef typename Kernel::source_type source_type;
    typedef typename Kernel::charge_type charge_type;
    typedef typename Kernel::target_type target_type;
    typedef typename Kernel::result_type result_type;

    static_assert(std::is_same<source_type,
                  typename std::iterator_traits<SourceIter>::value_type
                  >::value, "SourceIter::value_type != Kernel::source_type");
    static_assert(std::is_same<charge_type,
                  typename std::iterator_traits<ChargeIter>::value_type
                  >::value, "ChargeIter::value_type != Kernel::charge_type");
    static_assert(std::is_same<target_type,
                  typename std::iterator_traits<SourceIter>::value_type
                  >::value, "SourceIter::value_type != Kernel::target_type");
    static_assert(std::is_same<result_type,
                  typename std::iterator_traits<ResultIter>::value_type
                  >::value, "ResultIter::value_type != Kernel::result_type");

    for ( ; p1_first != p1_last; ++p1_first, ++c1_first, ++r1_first) {
      const source_type& p1 = *p1_first;
      const charge_type& c1 = *c1_first;
      result_type& r1       = *r1_first;

      SourceIter p2i = p2_first;
      ChargeIter c2i = c2_first;
      ResultIter r2i = r2_first;
      for ( ; p2i != p2_last; ++p2i, ++c2i, ++r2i)
        P2P::symm_eval(K, p1, c1, r1, *p2i, *c2i, *r2i);
    }
  }

  /** Symmetric diagonal block P2P using the evaluation operator
   * r_i += sum_j K(p_i, p_j) * c_j
   *
   * @pre source_type == target_type
   */
  template <typename Kernel,
            typename SourceIter, typename ChargeIter, typename ResultIter>
  inline static
  typename std::enable_if<!KernelTraits<Kernel>::has_vector_P2P_symm>::type
  block_eval(const Kernel& K,
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

    SourceIter pi = p_first;
    ChargeIter ci = c_first;
    ResultIter ri = r_first;
    // The first diagonal element
    *ri += K(*pi, *pi) * (*ci);

    for (++pi, ++ci, ++ri; pi != p_last; ++pi, ++ci, ++ri) {
      const source_type& p = *pi;
      const charge_type& c = *ci;
      result_type& r       = *ri;

      // The off-diagonal elements
      SourceIter pj = p_first;
      ChargeIter cj = c_first;
      ResultIter rj = r_first;
      for ( ; pj != pi; ++pj, ++cj, ++rj)
        P2P::symm_eval(K, p, c, r, *pj, *cj, *rj);

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
#if defined(FMMTL_DEBUG)
    std::cout << "P2P:"
              << "\n  " << source
              << "\n  " << target << std::endl;
#endif
    FMMTL_LOG("P2P 2box asymm");

    P2P::block_eval(c.kernel(),
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
#if defined(FMMTL_DEBUG)
    std::cout << "P2P:"
              << "\n  " << box1
              << "\n  " << box2 << std::endl;
    std::cout << "P2P:"
              << "\n  " << box2
              << "\n  " << box1 << std::endl;
#endif
    FMMTL_LOG("P2P 2box symm");

    P2P::block_eval(c.kernel(),
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
#if defined(FMMTL_DEBUG)
    std::cout << "P2P:"
              << "\n  " << box << std::endl;
#endif
    FMMTL_LOG("P2P 1box symm");

    P2P::block_eval(c.kernel(),
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

#include <unordered_map>

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
  typedef typename Context::source_box_type source_box_type;
  typedef typename Context::target_box_type target_box_type;

  std::vector<target_box_type> target_box_list;
  std::vector<std::vector<source_box_type>> source_boxes;
  unsigned box_pair_count;

  //! Box list for P2P interactions    TODO: could further compress these...
  //typedef std::pair<source_box_type, target_box_type> box_pair;
  //std::vector<box_pair> p2p_list;

  // For now, only use for GPU...
  P2P_Compressed<kernel_type>* p2p_compressed;

 public:
  P2P_Batch()
      : box_pair_count(0), p2p_compressed(nullptr) {
  }
  ~P2P_Batch() {
    delete p2p_compressed;
  }

  /** Insert a source-target box interaction to the interaction list */
  void insert(const source_box_type& s, const target_box_type& t) {
    if (source_boxes.size() <= t.index()) {
      source_boxes.resize(t.index() + 1);
      target_box_list.push_back(t);
    } else if (source_boxes[t.index()].empty()) {
      target_box_list.push_back(t);
    }

    source_boxes[t.index()].push_back(s);
    ++box_pair_count;

    //p2p_list.push_back(std::make_pair(s,t));
  }

  /** Compute all interations in the interaction list */
  void execute(Context& c) {
    FMMTL_LOG("P2P Batch");
#if FMMTL_NO_CUDA
    auto t_end = target_box_list.end();
#pragma omp parallel for
    for (auto ti = target_box_list.begin(); ti < t_end; ++ti) {
      target_box_type& tb = *ti;
      auto s_end = source_boxes[tb.index()].end();
      for (auto si = source_boxes[tb.index()].begin(); si != s_end; ++si) {
        P2P::eval(c, *si, tb, P2P::ONE_SIDED());
      }
    }
#else
    if (p2p_compressed == nullptr)
      p2p_compressed =
          P2P_Compressed<kernel_type>::make(c, target_box_list, source_boxes, box_pair_count);
    p2p_compressed->execute(c);
#endif
  }

  /*
  class P2P_Matrix
      : public EvaluatorBase<Context> {
    ublas::compressed_matrix<kernel_value_type> A;

   public:
    virtual void execute(Context& c) {
      // printf("EvalLocalSparse::execute(Context&)\n");

      typedef typename Context::charge_type charge_type;
      ublas::vector<charge_type> charges(c.source_tree().bodies());
      std::copy(c.charge_begin(), c.charge_end(), charges.begin());

      // Call the matvec
      typedef typename Context::result_type result_type;
      ublas::vector<result_type> results = ublas::prod(A, charges);

      // Accumulate results
      std::transform(results.begin(), results.end(),
                     c.result_begin(), c.result_begin(),
                     std::plus<result_type>());
    }
  };
  */


  /** Convert the interaction list to an interaction matrix
   * by evaluating all of the elements.
   */
  /*
  ublas::compressed_matrix<kernel_value_type> to_matrix(Context& bc) {
    auto first_source = bc.source_begin();
    auto first_target = bc.target_begin();

    // Interaction list for each target body
    unsigned rows = bc.target_tree().bodies();
    unsigned cols = 0;
    unsigned nnz = 0;
    std::vector<std::vector<unsigned>> csr(rows);

    for (const box_pair& b2b : p2p_list) {
      const source_box_type& box1 = b2b.first;
      const target_box_type& box2 = b2b.second;

      auto source_end = bc.source_end(box1);
      auto target_end = bc.target_end(box2);
      for (auto t = bc.target_begin(box2); t != target_end; ++t) {
        // Row
        unsigned i = t - first_target;
        std::vector<unsigned>& csri = csr[i];

        for (auto s = bc.source_begin(box1); s != source_end; ++s) {
          // Column
          unsigned j = s - first_source;

          //FMMTL_ASSERT(std::find(csri.begin(), csri.end(), j) == csri.end());
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
  */
};
