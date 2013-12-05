#pragma once
/** @file L2P.hpp
 * @brief Dispatch methods for the L2P stage
 *
 */

#include "fmmtl/Logger.hpp"
#include "fmmtl/meta/kernel_traits.hpp"
#include <type_traits>

struct L2P
{
  /** If no other L2P dispatcher matches */
  template <typename Expansion, typename... Args>
  inline static void eval(const Expansion&, Args...) {
    std::cerr << "Expansion does not have a correct L2P!\n";
    std::cerr << ExpansionTraits<Expansion>() << std::endl;
    exit(1);
  }

  /** L2P evaluation.
   * The Expansion provides a vector L2P accumulator.
   */
  template <typename Expansion, typename TargetIter, typename ResultIter>
  inline static
  typename std::enable_if<ExpansionTraits<Expansion>::has_vector_L2P>::type
  eval(const Expansion& K,
       const typename Expansion::local_type& L,
       const typename Expansion::point_type& center,
       TargetIter t_begin, TargetIter t_end, ResultIter r_begin) {
    K.L2P(L, center, t_begin, t_end, r_begin);
  }

  /** L2P evaluation.
   * The Expansion provides a scalar L2P accumulator. Use it for each target point.
   */
  template <typename Expansion, typename TargetIter, typename ResultIter>
  inline static
  typename std::enable_if<ExpansionTraits<Expansion>::has_scalar_L2P &
                          !ExpansionTraits<Expansion>::has_vector_L2P>::type
  eval(const Expansion& K,
       const typename Expansion::local_type& L,
       const typename Expansion::point_type& center,
       TargetIter t_begin, TargetIter t_end, ResultIter r_begin) {
    for ( ; t_begin != t_end; ++t_begin, ++r_begin)
      K.L2P(L, center, *t_begin, *r_begin);
  }

 public:

  /** Unwrap data from Context and dispatch to the L2P evaluator
   */
  template <typename Context>
  inline static void eval(Context& c,
                          const typename Context::target_box_type& tbox)
  {
#if defined(FMMTL_DEBUG)
    std::cout << "L2P:"
              << "\n  " << tbox << std::endl;
#endif
    FMMTL_LOG("L2P");

    L2P::eval(c.expansion(),
              c.local(tbox),
              tbox.center(),
              c.target_begin(tbox), c.target_end(tbox),
              c.result_begin(tbox));
  }
};
