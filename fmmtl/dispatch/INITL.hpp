#pragma once
/** @file INITL.hpp
 * @brief Dispatch methods for the initializing a local expansion
 *
 */

#include "fmmtl/Logger.hpp"
#include "fmmtl/meta/kernel_traits.hpp"
#include <type_traits>

struct INITL
{
  /** If no init_local dispatcher matches */
  template <typename Expansion, typename... Args>
  inline static void eval(const Expansion&,
                          typename Expansion::local_type& L,
                          Args...) {
    // Use default constructor
    L = typename Expansion::local_type();
  }

  template <typename Expansion>
  inline static
  typename std::enable_if<ExpansionTraits<Expansion>::has_init_local>::type
  eval(const Expansion& K,
       typename Expansion::local_type& L,
       const typename Expansion::point_type& extents,
       unsigned level) {
    K.init_local(L, extents, level);
  }

  /** Unwrap data from Context and dispatch to the INITL evaluator
   */
  template <typename Context>
  inline static void eval(Context& c,
                          const typename Context::target_box_type& tbox)
  {
#if defined(FMMTL_DEBUG)
    std::cout << "INITL:"
              << "\n  " << tbox << std::endl;
#endif
    FMMTL_LOG("INITL");

    INITL::eval(c.expansion(),
                c.local(tbox),
                tbox.extents(),
                tbox.level());
  }
};
