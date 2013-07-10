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

#include "executor/P2P_Compressed.hpp"

#if defined(FMMTL_KERNEL)
#if defined(FMMTL_NO_CUDA)
#include "executor/P2P_Compressed.cpp"
#endif
#if defined(FMMTL_USE_CUDA) && defined(__CUDACC__)
#include "executor/P2P_Compressed.cu"
#endif
#else
#endif

template <class Kernel>
struct FMM_Kernel {
  typedef Kernel kernel_type;
  // TODO: Do the template instantiations for linking automatically...?
  // TODO: Vectorized and/or GPU evaluations for Kernels?
};

// Template instantiations for external compilation and linking
// TODO: Remove?
#define FMMTL_KERNEL_EXTRAS(kernel) template class P2P_Compressed<kernel>


template <class Expansion>
struct FMM_Expansion {
  typedef Expansion expansion_type;

  // TODO: Add syntactic sugar for Expansions
  // eg fmm_matrix<expansion_type> operator()(std::vector<source_type>,
  //                                        std::vector<target_type>);
};
