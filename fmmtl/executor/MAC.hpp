#pragma once
/** @file MAC.hpp
 * @brief Dispatch methods for the Multipole Acceptance Criteria
 */

#include "fmmtl/Logger.hpp"
#include "fmmtl/KernelTraits.hpp"
#include <type_traits>

class MAC
{
  /** If no other MAC dispatcher matches */
  template <typename Expansion, typename... Args>
  inline static void eval(const Expansion&, Args...) {
    // TODO
  }

  template <typename Expansion>
  inline static
  typename std::enable_if<ExpansionTraits<Expansion>::has_M2L>::type
  eval(const Expansion& K,
       const typename Expansion::multipole_type& source,
       typename Expansion::local_type& target,
       const typename Expansion::point_type& translation) {
    // TODO
  }

 public:

  /** Unwrap data from Context and dispatch to the MAC evaluator
   */
  template <typename Context>
  inline static void eval(Context& c,
                          const typename Context::source_box_type& sbox,
                          const typename Context::target_box_type& tbox)
  {
#if defined(FMMTL_DEBUG)
    std::cout << "MAC:"
              << "\n  " << sbox
              << "\n  " << tbox << std::endl;
#endif
    FMMTL_LOG("MAC");

    MAC::eval(c.expansion(),
              c.multipole(sbox),
              c.local(tbox),
              tbox.center() - sbox.center());
  }
};
