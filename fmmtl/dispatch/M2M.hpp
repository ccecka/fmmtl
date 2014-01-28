#pragma once
/** @file M2M.hpp
 * @brief Dispatch methods for the M2M stage
 *
 */

#include "fmmtl/Logger.hpp"
#include "fmmtl/meta/kernel_traits.hpp"
#include <type_traits>

class M2M
{
  /** If no other M2M dispatcher matches */
  template <typename Expansion, typename... Args>
  inline static void eval(const Expansion&, Args...) {
    std::cerr << "Expansion does not have a correct M2M!\n";
    std::cerr << ExpansionTraits<Expansion>() << std::endl;
    exit(1);
  }

  template <typename Expansion>
  inline static
  typename std::enable_if<ExpansionTraits<Expansion>::has_M2M>::type
  eval(const Expansion& K,
       const typename Expansion::multipole_type& source,
       typename Expansion::multipole_type& target,
       const typename Expansion::point_type& translation) {
    K.M2M(source, target, translation);
  }

 public:

  /** Unwrap data from Context and dispatch to the M2M evaluator
   */
  template <typename Context>
  inline static void eval(Context& c,
                          const typename Context::source_box_type& sbox,
                          const typename Context::source_box_type& tbox)
  {
#if defined(FMMTL_DEBUG)
    std::cout << "M2M:"
              << "\n  " << sbox
              << "\n  " << tbox << std::endl;
#endif
    FMMTL_LOG("M2M");

    M2M::eval(c.expansion(),
              c.multipole(sbox),
              c.multipole(tbox),
              tbox.center() - sbox.center());
  }
};
