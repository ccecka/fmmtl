#pragma once

#if defined(__CUDACC__)
#  define FMMTL_WITH_CUDA
#endif

// FMMTL_INLINE
#if defined(__CUDACC__)        // Compiling with nvcc
#  define FMMTL_INLINE __host__ __device__ inline
#else                          // Not compiling with nvcc
#  define FMMTL_INLINE inline
#endif

#include <thrust/version.h>
#if (THRUST_VERSION < 100700)
#  error Need Thrust v1.7. Please upgrade to CUDA 5.5.
#endif

// Enable performance options in NDEBUG mode
#if defined(FMMTL_NDEBUG) || defined(NDEBUG)
#  define FMMTL_DISABLE_ASSERTS
#  define FMMTL_DISABLE_CUDA_CHECKS
#endif

// Disable performance options in DEBUG mode
#if defined(FMMTL_DEBUG)
#  undef FMMTL_DISABLE_ASSERTS
#  undef FMMTL_DISABLE_CUDA_CHECKS
#endif

// FMMTL_ASSERT
#undef FMMTL_ASSERT
#if defined(FMMTL_DISABLE_ASSERTS)
#  define FMMTL_ASSERT(expr) ((void)0)
#else
#  include <cassert>
#  define FMMTL_ASSERT(expr) assert(expr)
#endif

// FMMTL_CUDA_CHECK
#undef FMMTL_CUDA_CHECK
#if !defined(FMMTL_WITH_CUDA) || defined(FMMTL_DISABLE_CUDA_CHECKS)
#  define FMMTL_CUDA_CHECK ((void)0)
#else
#  include <iostream>
#  include <thrust/system/cuda/error.h>
inline void cuda_check(const char* file, int line) {
  cudaError_t error = cudaDeviceSynchronize();
  if (error != cudaSuccess) {
    std::cerr << "CUDA assert:"
              << " " << cudaGetErrorString(error)
              << " " << file
              << ":" << line
              << std::endl;;
    exit(error);
  }
}
#  define FMMTL_CUDA_CHECK cuda_check(__FILE__, __LINE__)
#endif

// Load OpenMP support if available
#if defined(_OPENMP)
#  include <omp.h>
#else
#  warning Compiler does not support OpenMP
typedef int omp_int_t;
inline omp_int_t omp_get_thread_num()  { return 0; }
inline omp_int_t omp_get_max_threads() { return 1; }
#endif
