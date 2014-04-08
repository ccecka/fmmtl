#pragma once
/** @file M2P.hpp
 * @brief Dispatch methods for the M2P stage
 *
 */

#include "fmmtl/util/Logger.hpp"
#include "fmmtl/meta/kernel_traits.hpp"
#include <type_traits>

/** Default behavior gives a warning -- using non-existent method */
template <bool has_m2p>
struct M2P_Helper {
  inline static void apply(...) {
    std::cerr << "WARNING: Expansion does not have a correct M2P!\n";
  }
  inline static void eval(...) {
    std::cerr << "WARNING: Expansion does not have a correct M2P!\n";
  }
};

/** Expansion has an M2P method (scalar/vector) to dispatch to */
template <>
struct M2P_Helper<true> {
  /** The Expansion provides a vector M2P accumulator. */
  template <typename Expansion, typename TargetIter, typename ResultIter>
  inline static
  typename std::enable_if<ExpansionTraits<Expansion>::has_vector_M2P>::type
  apply(const Expansion& K,
        const typename Expansion::multipole_type& M,
        const typename Expansion::point_type& center,
        TargetIter t_begin, TargetIter t_end, ResultIter r_begin) {
    K.M2P(M, center, t_begin, t_end, r_begin);
  }

  /** The Expansion provides a scalar M2P accumulator. */
  template <typename Expansion>
  inline static
  typename std::enable_if<ExpansionTraits<Expansion>::has_scalar_M2P &
                          !ExpansionTraits<Expansion>::has_vector_M2P>::type
  apply(const Expansion& K,
        const typename Expansion::multipole_type& M,
        const typename Expansion::point_type& center,
        const typename Expansion::target_type& target,
              typename Expansion::result_type& result) {
    K.M2P(M, center, target, result);
  }

  /** The Expansion provides a scalar M2P accumulator. */
  template <typename Expansion, typename TargetIter, typename ResultIter>
  inline static
  typename std::enable_if<ExpansionTraits<Expansion>::has_scalar_M2P &
                          !ExpansionTraits<Expansion>::has_vector_M2P>::type
  apply(const Expansion& K,
        const typename Expansion::multipole_type& M,
        const typename Expansion::point_type& center,
        TargetIter t_begin, TargetIter t_end, ResultIter r_begin) {
    for ( ; t_begin != t_end; ++t_begin, ++r_begin)
      apply(K, M, center, *t_begin, *r_begin);
  }

  /** Unpack from Context and apply */
  template <typename Context>
  inline static void eval(Context& c,
                          const typename Context::source_box_type& sbox,
                          const typename Context::target_box_type& tbox) {
    apply(c.expansion(),
          c.multipole(sbox), sbox.center(),
          c.target_begin(tbox), c.target_end(tbox), c.result_begin(tbox));
  }
};

/** Public M2P dispatcher */
struct M2P {
  /** Forward to M2P_Helper::apply */
  template <typename Expansion>
  inline static void apply(const Expansion& K,
                           const typename Expansion::multipole_type& M,
                           const typename Expansion::point_type& center,
                           const typename Expansion::target_type& target,
                           typename Expansion::result_type& result) {
    typedef M2P_Helper<ExpansionTraits<Expansion>::has_M2P> M2P_H;
    M2P_H::apply(K, M, center, target, result);
  }

  /** Forward to M2P_Helper::eval */
  template <typename Context>
  inline static void eval(Context& c,
                          const typename Context::source_box_type& sbox,
                          const typename Context::target_box_type& tbox)
  {
#if defined(FMMTL_DEBUG)
    std::cout << "M2P:"
              << "\n  " << sbox
              << "\n  " << tbox << std::endl;
#endif
    FMMTL_LOG("M2P");

    typedef ExpansionTraits<typename Context::expansion_type> expansion_traits;
    typedef M2P_Helper<expansion_traits::has_M2P> M2P_H;
    M2P_H::eval(c, sbox, tbox);
  }
};
