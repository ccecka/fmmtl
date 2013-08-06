#pragma once


#if defined(FMMTL_NO_CUDA)    // Disable CUDA/Thrust acceleration

# if defined(__CUDACC__)
#  error Compiling with nvcc and NO_CUDA flag!
# endif

# undef FMMTL_USE_THRUST

#else                         // Enable CUDA/Thrust acceleration

# if defined(FMMTL_KERNEL) & !defined(__CUDACC__)
#  error Should be compiling .kern with CUDA!
# endif

# define FMMTL_USE_THRUST
# include <thrust/version.h>
# if (THRUST_VERSION < 100700)
#  error Need Thrust 1.7. Please upgrade to CUDA 5.5.
# endif
# include <thrust/detail/config.h>

# if defined(__CUDACC__)        // Compiling with nvcc
#  define FMMTL_INLINE inline __host__ __device__   // Provide CUDA methods
#  define FMMTL_KERNEL                              // nvcc only for .kern now
# endif
#endif  // end FMMTL_NO_CUDA


#if !defined(FMMTL_INLINE)
# define FMMTL_INLINE inline
#endif

// Enable performance options in NDEBUG mode
#if defined(FMMTL_NDEBUG)
# define FMMTL_DISABLE_ASSERTS
#endif

// Disable performance options in DEBUG mode
#if defined(FMMTL_DEBUG)
# undef FMMTL_DISABLE_ASSERTS
#endif

// FMMTL assert
#undef FMMTL_ASSERT
#if defined(FMMTL_DISABLE_ASSERTS)
# define FMMTL_ASSERT(expr) ((void)0)
#else
# include <cassert>
# define FMMTL_ASSERT(expr) assert(expr)
#endif
