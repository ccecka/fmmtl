#pragma once
/** @file FMM_Kernel
 *
 * CRTP base classes for Kernels and Expansions.
 * This allows the future extension of Kernels and Expansions with methods that
 * need not be implemented within the Kernel/Expansion.
 * Additionally, the traits system could be incorporated here to give the
 * library improved access and insight to these classes.
 *
 * At the moment, these may not be fully necessary and are subject to change.
 */


#include "config.hpp"

#if defined(FMMTL_KERNEL)   // If compiling the .kern, include implementations
#if defined(FMMTL_NO_CUDA)                          // If not using CUDA
#include "executor/P2P_Compressed.cpp"
#endif
#if defined(FMMTL_USE_CUDA) && defined(__CUDACC__)  // If compiling with nvcc
#include "executor/P2P_Compressed.cu"
#endif
#endif

#if !defined(FMMTL_KERNEL)  // If !compiling the .kern, include headers and C++11
#include "executor/P2P_Compressed.hpp"

#include "KernelMatrix.hpp"
#endif


template <class Kernel>
struct FMM_Kernel {
  typedef Kernel kernel_type;

  Kernel& kernel() {
    return static_cast<Kernel&>(*this);
  }
  const Kernel& kernel() const {
    return static_cast<const Kernel&>(*this);
  }

  // TODO: Do the template instantiations for linking automatically...?
  // TODO: Vectorized and/or GPU evaluations for Kernels?
};

// Template instantiations for external compilation and linking
// TODO: Remove?
#define FMMTL_KERNEL_EXTRAS(kernel) template class P2P_Compressed<kernel>


template <class Kernel, class Expansion>
class FMM_Expansion
    : public Kernel {
 public:
  typedef typename Kernel::target_type target_type;
  typedef typename Kernel::source_type source_type;

  Expansion& expansion() {
    return static_cast<Expansion&>(*this);
  }
  const Expansion& expansion() const {
    return static_cast<const Expansion&>(*this);
  }
  Kernel& kernel() {
    return static_cast<Kernel&>(*this);
  }
  const Kernel& kernel() const {
    return static_cast<const Kernel&>(*this);
  }

  class KernelMatrixProxy {
    const Expansion& E_;
    const std::vector<target_type>& t_;
    const std::vector<source_type>& s_;

   public:
    KernelMatrixProxy(const Expansion& E,
                      const std::vector<target_type>& t,
                      const std::vector<source_type>& s)
        : E_(E), t_(t), s_(s)  {
    }

#if !defined(FMMTL_KERNEL)
    operator fmm_matrix<Expansion>() {
      if (&s_ == &t_)
        return make_fmm_matrix(E_, s_, t_);
      else
        return make_fmm_matrix(E_, s_);
    }
#endif
  };

  KernelMatrixProxy operator()(const std::vector<target_type>& t,
                               const std::vector<source_type>& s) const {
    return KernelMatrixProxy(expansion(), t, s);
  }
};
