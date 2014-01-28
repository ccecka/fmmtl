#pragma once
/** @file M2P.hpp
 * @brief Dispatch methods for the M2P stage
 *
 */

#include "fmmtl/Logger.hpp"
#include "fmmtl/meta/kernel_traits.hpp"
#include <type_traits>

class P2L
{
  /** If no other P2L dispatcher matches */
  template <typename Expansion, typename... Args>
  inline static void eval(const Expansion&, Args...) {
    std::cerr << "Expansion does not have a correct P2L!\n";
    std::cerr << ExpansionTraits<Expansion>() << std::endl;
    exit(1);
  }

  /** P2L evaluation.
   * The Expansion provides a vector P2L accumulator.
   */
  template <typename Expansion, typename SourceIter, typename ChargeIter>
  inline static
  typename std::enable_if<ExpansionTraits<Expansion>::has_vector_P2L>::type
  eval(const Expansion& K,
       SourceIter s_begin, SourceIter s_end,
       ChargeIter c_begin,
       const typename Expansion::point_type& center,
       typename Expansion::local_type& L) {
    K.P2L(s_begin, s_end, c_begin, center, L);
  }

  /** P2L evaluation.
   * The Expansion provides a scalar P2L accumulator. Use it for each source.
   */
  template <typename Expansion, typename SourceIter, typename ChargeIter>
  inline static
  typename std::enable_if<ExpansionTraits<Expansion>::has_scalar_P2L &
                          !ExpansionTraits<Expansion>::has_vector_P2L>::type
  eval(const Expansion& K,
       SourceIter s_begin, SourceIter s_end,
       ChargeIter c_begin,
       const typename Expansion::point_type& center,
       typename Expansion::local_type& L) {
    for ( ; s_begin != s_end; ++s_begin, ++c_begin)
      K.P2L(*s_begin, *c_begin, center, L);
  }

 public:

  /** Unwrap the data from BoxContext and dispatch to the P2L evaluator
   */
  template <typename Context>
  inline static void eval(Context& c,
                          const typename Context::source_body_iterator sfirst,
                          const typename Context::source_body_iterator slast,
                          const typename Context::target_box_type& tbox) {
#if defined(FMMTL_DEBUG)
    std::cout << "P2L:"
              << "\n  Bodies [" << sfirst << ", " << slast << ")"
              << "\n  " << tbox << std::endl;
#endif
    FMMTL_LOG("P2L vec");

    P2L::eval(c.expansion(),
              c.source(sfirst), c.source(slast),
              c.charge(sfirst),
              tbox.center(),
              c.local(tbox));
  }

  /** Unwrap the data from BoxContext and dispatch to the P2L evaluator
   */
  template <typename Context>
  inline static void eval(Context& c,
                          const typename Context::source_box_type& source,
                          const typename Context::target_box_type& target) {
#if defined(FMMTL_DEBUG)
    std::cout << "P2L:"
              << "\n  " << source
              << "\n  " << target << std::endl;
#endif
    FMMTL_LOG("P2L box");

    P2L::eval(c.expansion(),
              c.source_begin(source), c.source_begin(source),
              c.charge_begin(source),
              target.center(),
              c.local(target));
  }
};
