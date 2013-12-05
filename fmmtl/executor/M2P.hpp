#pragma once
/** @file M2P.hpp
 * @brief Dispatch methods for the M2P stage
 *
 */

#include "fmmtl/Logger.hpp"
#include "fmmtl/meta/kernel_traits.hpp"
#include <type_traits>

class M2P
{
  /** If no other M2P dispatcher matches */
  template <typename Expansion, typename... Args>
  inline static void eval(const Expansion&, Args...) {
    std::cerr << "Expansion does not have a correct M2P!\n";
    std::cerr << ExpansionTraits<Expansion>() << std::endl;
    exit(1);
  }

  /** M2P evaluation.
   * The Expansion provides a vector M2P accumulator.
   */
  template <typename Expansion, typename TargetIter, typename ResultIter>
  inline static
  typename std::enable_if<ExpansionTraits<Expansion>::has_vector_M2P>::type
  eval(const Expansion& K,
       const typename Expansion::multipole_type& M,
       const typename Expansion::point_type& center,
       TargetIter t_begin, TargetIter t_end,
       ResultIter r_begin) {
    K.M2P(M, center, t_begin, t_end, r_begin);
  }

  /** M2P evaluation.
   * The Expansion provides a scalar M2P accumulator. Use it for each target point.
   */
  template <typename Expansion, typename TargetIter, typename ResultIter>
  inline static
  typename std::enable_if<ExpansionTraits<Expansion>::has_scalar_M2P &
                          !ExpansionTraits<Expansion>::has_vector_M2P>::type
  eval(const Expansion& K,
       const typename Expansion::multipole_type& M,
       const typename Expansion::point_type& center,
       TargetIter t_begin, TargetIter t_end,
       ResultIter r_begin) {
    for ( ; t_begin != t_end; ++t_begin, ++r_begin)
      K.M2P(M, center, *t_begin, *r_begin);
  }

 public:

  /** Unwrap data from Context and dispatch to the M2P evaluator
   */
  template <typename Context>
  inline static void eval(Context& c,
                          const typename Context::source_box_type& sbox,
                          const typename Context::target_body_iterator tfirst,
                          const typename Context::target_body_iterator tlast) {
#if defined(FMMTL_DEBUG)
    std::cout << "M2P:"
              << "\n  " << sbox
              << "\n  Bodies [" << tfirst << ", " << tlast << ")" << std::endl;
#endif
    FMMTL_LOG("M2P vec");

    M2P::eval(c.expansion(),
              c.multipole(sbox),
              sbox.center(),
              c.target(tfirst), c.target(tlast),
              c.result(tfirst));
  }

  /** Unwrap data from Context and dispatch to the M2P evaluator
   */
  template <typename Context>
  inline static void eval(Context& c,
                          const typename Context::source_box_type& source,
                          const typename Context::target_box_type& target) {
#if defined(FMMTL_DEBUG)
    std::cout << "M2P:"
              << "\n  " << source
              << "\n  " << target << std::endl;
#endif
    FMMTL_LOG("M2P box");

    M2P::eval(c.expansion(),
              c.multipole(source),
              source.center(),
              c.target_begin(target), c.target_end(target),
              c.result_begin(target));
  }
};
