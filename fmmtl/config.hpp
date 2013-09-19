#pragma once


#if defined(FMMTL_NO_CUDA)    // Disable CUDA/Thrust acceleration

# if defined(__CUDACC__)
#  error Compiling with nvcc and NO_CUDA flag!
# endif

#else                         // Enable CUDA/Thrust acceleration

# if defined(FMMTL_KERNEL) & !defined(__CUDACC__)
#  error Should be compiling .kern with CUDA!
# endif

# include <thrust/version.h>
# if (THRUST_VERSION < 100700)
#  error Need Thrust v1.7. Please upgrade to CUDA 5.5.
# endif
# include <thrust/detail/config.h>
# include <cublas.h>

#endif  // end FMMTL_NO_CUDA

// FMMTL_INLINE
#if defined(__CUDACC__)        // Compiling with nvcc
# define FMMTL_INLINE inline __host__ __device__   // Provide CUDA methods
#else
# define FMMTL_INLINE inline
#endif

// Enable performance options in NDEBUG mode
#if defined(FMMTL_NDEBUG)
# define FMMTL_DISABLE_ASSERTS
# define FMMTL_DISABLE_CUDA_CHECKS
#endif

// Disable performance options in DEBUG mode
#if defined(FMMTL_DEBUG)
# undef FMMTL_DISABLE_ASSERTS
# undef FMMTL_DISABLE_CUDA_CHECKS
#endif

// FMMTL_ASSERT
#undef FMMTL_ASSERT
#if defined(FMMTL_DISABLE_ASSERTS)
# define FMMTL_ASSERT(expr) ((void)0)
#else
# include <cassert>
# define FMMTL_ASSERT(expr) assert(expr)
#endif

// FMMTL_CUDA_CHECK
#undef FMMTL_CUDA_CHECK
#if defined(FMMTL_DISABLE_CUDA_CHECKS) || defined(FMMTL_NO_CUDA)
# define FMMTL_CUDA_CHECK
#else
#include <cstdlib>
#include <cstdio>
inline void cuda_check(char* file, int line) {
  cudaError_t code = cudaDeviceSynchronize();
  if (code != cudaSuccess) {
    fprintf(stderr,"CUDA assert: %s %s %d\n",cudaGetErrorString(code),file,line);
    exit(code);
  }
}
# define FMMTL_CUDA_CHECK cuda_check(__FILE__, __LINE__)
#endif
