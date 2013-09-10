#pragma once
/** @file Expansion
 *
 * CRTP base class for Expansions.
 * This allows the future extension of Expansions with methods that
 * need not be implemented within the Expansion.
 * Additionally, the traits system could be incorporated here to give the
 * library improved access and insight to these classes.
 *
 * At the moment, these may not be fully necessary and are subject to change.
 */

#include "config.hpp"

// Protect against Expansions defined within .kern files
#if !defined(FMMTL_KERNEL)
# include "KernelMatrix.hpp"
#endif

namespace fmmtl
{

template <class Kernel, class DerivedExpansion>
struct Expansion
    : public Kernel {
  // The type of kernel this expansion is representing
  typedef Kernel            kernel_type;
  // The type of the derived expansion implementation
  typedef DerivedExpansion  expansion_type;

  // Default constructor
  Expansion() {}

  // Forward the construction on a Kernel to the copy constructor
  Expansion(const Kernel& K)
      : Kernel(K) {}

  expansion_type& expansion() {
    return static_cast<expansion_type&>(*this);
  }
  const expansion_type& expansion() const {
    return static_cast<const expansion_type&>(*this);
  }
  kernel_type& kernel() {
    return static_cast<kernel_type&>(*this);
  }
  const kernel_type& kernel() const {
    return static_cast<const kernel_type&>(*this);
  }

#if !defined(FMMTL_KERNEL)
  typedef typename kernel_type::target_type target_type;
  typedef typename kernel_type::source_type source_type;

  // Sugar that can be improved
  class KernelMatrixProxy {
    const expansion_type& E_;
    const std::vector<target_type>& t_;
    const std::vector<source_type>& s_;

   public:
    KernelMatrixProxy(const expansion_type& E,
                      const std::vector<target_type>& t,
                      const std::vector<source_type>& s)
        : E_(E), t_(t), s_(s)  {
    }

    operator kernel_matrix<expansion_type>() {
      if (&s_ == &t_)
        return make_kernel_matrix(E_, s_, t_);
      else
        return make_kernel_matrix(E_, s_);
    }
  };

  KernelMatrixProxy operator()(const std::vector<target_type>& t,
                               const std::vector<source_type>& s) const {
    return KernelMatrixProxy(expansion(), t, s);
  }
#endif
};

} // end namespace fmmtl
