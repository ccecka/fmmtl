#pragma once
/** @file INITL.hpp
 * @brief Dispatch methods for the initializing a local expansion
 *
 */

#include "fmmtl/KernelTraits.hpp"
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

  template <typename Context>
  inline static void eval(Context& c,
                          const typename Context::target_box_type& tbox)
  {
#ifdef DEBUG
    std::cout << "INITL: " << tbox << std::endl;
#endif

    INITL::eval(c.expansion(), c.local(tbox), tbox.extents(), tbox.level());
  }
};
