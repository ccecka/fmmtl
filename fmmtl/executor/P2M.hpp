#pragma once
/** @file P2M.hpp
 * @brief Dispatch methods for the P2M stage
 *
 */

#include "fmmtl/Logger.hpp"
#include "fmmtl/KernelTraits.hpp"
#include <type_traits>

class P2M
{
  /** If no other P2M dispatcher matches */
  template <typename Expansion, typename... Args>
  inline static void eval(const Expansion&, Args...) {
    std::cerr << "Expansion does not have a correct P2M!\n";
    std::cerr << ExpansionTraits<Expansion>() << std::endl;
    exit(1);
  }

  /** P2M evaluation.
   * The Expansion provides a vector P2M accumulator.
   */
  template <typename Expansion, typename SourceIter, typename ChargeIter>
  inline static
  typename std::enable_if<ExpansionTraits<Expansion>::has_vector_P2M>::type
  eval(const Expansion& K,
       SourceIter s_begin, SourceIter s_end,
       ChargeIter c_begin,
       const typename Expansion::point_type& center,
       typename Expansion::multipole_type& M) {
    K.P2M(s_begin, s_end, c_begin, center, M);
  }

  /** P2M evaluation.
   * The Expansion provides a scalar P2M accumulator. Use it for each source point.
   */
  template <typename Expansion, typename SourceIter, typename ChargeIter>
  inline static
  typename std::enable_if<ExpansionTraits<Expansion>::has_scalar_P2M &
                          !ExpansionTraits<Expansion>::has_vector_P2M>::type
  eval(const Expansion& K,
       SourceIter s_begin, SourceIter s_end,
       ChargeIter c_begin,
       const typename Expansion::point_type& center,
       typename Expansion::multipole_type& M) {
    for ( ; s_begin != s_end; ++s_begin, ++c_begin)
      K.P2M(*s_begin, *c_begin, center, M);
  }

 public:

  /** Unwrap the data from the Context and dispatch to P2M evaluator
   */
  template <typename Context>
  inline static void eval(Context& c,
                          const typename Context::source_box_type& sbox)
  {
#if defined(FMMTL_DEBUG)
    std::cout << "P2M:"
              << "\n  " << sbox << std::endl;
#endif
    FMMTL_LOG("P2M");

    P2M::eval(c.expansion(),
              c.source_begin(sbox), c.source_end(sbox),
              c.charge_begin(sbox),
              sbox.center(),
              c.multipole(sbox));
  }
};
