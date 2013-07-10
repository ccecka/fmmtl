#pragma once
/** @file M2P.hpp
 * @brief Dispatch methods for the M2P stage
 *
 */

#include "KernelTraits.hpp"

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
  typename std::enable_if<ExpansionTraits<Expansion>::has_M2P &
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

  /** Unwrap the data from BoxContext and dispatch to the M2P evaluator
   */
  template <typename Context>
  inline static void eval(Context& c,
                          const typename Context::source_box_type& source,
                          const typename Context::target_body_iterator tfirst,
                          const typename Context::target_body_iterator tlast) {
#ifdef DEBUG
    printf("M2P: %d to [%d,%d)\n", source.index(), tfirst.index(), tlast.index());
#endif

    M2P::eval(c.expansion(),
              c.multipole(source),
              source.center(),
              c.target(tfirst), c.target(tlast),
              c.result(tfirst));
  }

  /** Unwrap the data from BoxContext and dispatch to the M2P evaluator
   */
  template <typename Context>
  inline static void eval(Context& c,
                          const typename Context::source_box_type& source,
                          const typename Context::target_box_type& target) {
#ifdef DEBUG
    printf("M2P: %d to %d\n", source.index(), target.index());
#endif

    M2P::eval(c.expansion(),
              c.multipole(source),
              source.center(),
              c.target_begin(target), c.target_end(target),
              c.result_begin(target));
  }
};
