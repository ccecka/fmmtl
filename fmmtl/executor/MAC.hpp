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
  inline static bool eval(const Expansion&, Args...) {
    std::cerr << "Expansion does not have a correct dynamic MAC!\n";
    std::cerr << ExpansionTraits<Expansion>() << std::endl;
    exit(1);
  }

  template <typename Expansion>
  inline static
  typename std::enable_if<ExpansionTraits<Expansion>::has_dynamic_MAC,bool>::type
  eval(const Expansion& K,
       const typename Expansion::multipole_type& M,
       const typename Expansion::local_type& L) {
    return K.MAC(M, L);
  }

 public:

  /** Unwrap data from Context and dispatch to the MAC evaluator
   */
  template <typename Context>
  inline static bool eval(Context& c,
                          const typename Context::source_box_type& sbox,
                          const typename Context::target_box_type& tbox)
  {
#if defined(FMMTL_DEBUG)
    std::cout << "MAC:"
              << "\n  " << sbox
              << "\n  " << tbox << std::endl;
#endif
    FMMTL_LOG("MAC");

    bool dyn_mac = true;
    // Avoid asking for the multipole/local if we don't need them
    if (ExpansionTraits<typename Context::expansion_type>::has_dynamic_MAC) {
      dyn_mac = MAC::eval(c.expansion(),
                          c.multipole(sbox),
                          c.local(tbox));
    }

    return dyn_mac && c.mac(sbox,tbox);
  }
};
