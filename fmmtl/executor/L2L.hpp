#pragma once
/** @file L2L.hpp
 * @brief Dispatch methods for the L2L stage
 *
 */

#include "fmmtl/Logger.hpp"
#include "fmmtl/KernelTraits.hpp"
#include <type_traits>

struct L2L
{
  /** If no other L2L dispatcher matches */
  template <typename Expansion, typename... Args>
  inline static void eval(const Expansion&, Args...) {
    std::cerr << "Expansion does not have a correct L2L!\n";
    std::cerr << ExpansionTraits<Expansion>() << std::endl;
    exit(1);
  }

  template <typename Expansion>
  inline static
  typename std::enable_if<ExpansionTraits<Expansion>::has_L2L>::type
  eval(const Expansion& K,
       const typename Expansion::local_type& source,
       typename Expansion::local_type& target,
       const typename Expansion::point_type& translation) {
    K.L2L(source, target, translation);
  }

 public:

  /** Unwrap data from Context and dispatch to the L2L evaluator
   */
  template <typename Context>
  inline static void eval(Context& c,
                          const typename Context::target_box_type& sbox,
                          const typename Context::target_box_type& tbox)
  {
#ifdef DEBUG
    std::cout << "L2L:"
              << "\n  " << sbox
              << "\n  " << tbox << std::endl;
#endif
    FMMTL_LOG("L2L");

    L2L::eval(c.expansion(),
              c.local(sbox),
              c.local(tbox),
              tbox.center() - sbox.center());
  }
};

