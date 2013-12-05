#pragma once
/** @file INITM.hpp
 * @brief Dispatch methods for the initializing a multipole expansion
 *
 */

#include "fmmtl/Logger.hpp"
#include "fmmtl/meta/kernel_traits.hpp"
#include <type_traits>

struct INITM
{
  /** If no init_multipole dispatcher matches */
  template <typename Expansion, typename... Args>
  inline static void eval(const Expansion&,
                          typename Expansion::multipole_type& M,
                          Args...) {
    // Use default constructor
    M = typename Expansion::multipole_type();
  }

  template <typename Expansion>
  inline static
  typename std::enable_if<ExpansionTraits<Expansion>::has_init_multipole>::type
  eval(const Expansion& K,
       typename Expansion::multipole_type& M,
       const typename Expansion::point_type& extents,
       unsigned level) {
    K.init_multipole(M, extents, level);
  }

  /** Unwrap data from Context and dispatch to the INITM evaluator
   */
  template <typename Context>
  inline static void eval(Context& c,
                          const typename Context::source_box_type& sbox)
  {
#if defined(FMMTL_DEBUG)
    std::cout << "INITM:"
              << "\n  " << sbox << std::endl;
#endif
    FMMTL_LOG("INITM");

    INITM::eval(c.expansion(),
                c.multipole(sbox),
                sbox.extents(),
                sbox.level());
  }
};
