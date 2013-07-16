#pragma once
/** @file M2L.hpp
 * @brief Dispatch methods for the M2L stage
 *
 */

#include "fmmtl/KernelTraits.hpp"
#include <type_traits>

class M2L
{
  /** If no other M2L dispatcher matches */
  template <typename Expansion, typename... Args>
  inline static void eval(const Expansion&, Args...) {
    std::cerr << "Expansion does not have a correct M2L!\n";
    std::cerr << ExpansionTraits<Expansion>() << std::endl;
    exit(1);
  }

  template <typename Expansion>
  inline static
  typename std::enable_if<ExpansionTraits<Expansion>::has_M2L>::type
  eval(const Expansion& K,
       const typename Expansion::multipole_type& source,
       typename Expansion::local_type& target,
       const typename Expansion::point_type& translation) {
    K.M2L(source, target, translation);
  }

 public:

  template <typename Context>
  inline static void eval(Context& c,
                          const typename Context::source_box_type& source,
                          const typename Context::target_box_type& target)
  {
#ifdef DEBUG
    std::cout << "M2L:\n  " << source << "\n  " << target << std::endl;
#endif

    M2L::eval(c.expansion(),
              c.multipole(source),
              c.local(target),
              target.center() - source.center());
  }
};
