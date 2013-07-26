#pragma once

// Disable CUDA/Thrust acceleration
#if defined(FMMTL_NO_CUDA)

# undef FMMTL_USE_THRUST
# if defined(__CUDACC__)
#  error Compiling with nvcc and NO_CUDA flag!
# endif

// Enable CUDA/Thrust acceleration
#else

# define FMMTL_USE_THRUST
# if defined(__CUDACC__)        // Compiling with nvcc
#  define FMMTL_INLINE inline __host__ __device__
# endif

#endif


#if !defined(FMMTL_INLINE)
# define FMMTL_INLINE inline
#endif


// Enable performance options in NDEBUG mode
#if defined(FMMTL_NDEBUG)
# define FMMTL_DISABLE_ASSERTS
#endif


// Disable performance options in DEBUG mode
#if defined(FMMTL_DEBUG)
#define FMMTL_INLINE
#endif


#undef FMMTL_ASSERT

#if defined(FMMTL_DISABLE_ASSERTS)
# define FMMTL_ASSERT(expr) ((void)0)
#else
# include <cassert>
# define FMMTL_ASSERT(expr) assert(expr)
#endif
