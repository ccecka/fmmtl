#pragma once

// Disable CUDA/Thrust acceleration
#if defined(FMMTL_NO_CUDA)

// Enable CUDA/Thrust acceleration
#else

#define FMMTL_USE_THRUST 1
#define FMMTL_USE_CUDA 1
#ifdef __CUDACC__        // Compiling with nvcc
#define FMMTL_INLINE inline __host__ __device__
#endif

#endif


#if defined(FMMTL_NO_CUDA) && defined(FMMTL_USE_CUDA)
#error Undefined CUDA usage
#endif


// Enable performance options in RELEASE mode
#if defined(NDEBUG) || defined(FMMTL_NDEBUG)

#ifndef FMMTL_INLINE
#define FMMTL_INLINE inline
#endif
#ifndef FMMTL_CHECK
#define FMMTL_CHECK 0
#endif

// Disable performance options in DEBUG mode
#else

#ifndef FMMTL_INLINE
#define FMMTL_INLINE
#endif
#ifndef FMMTL_CHECK
#define FMMTL_CHECK 1
#endif

#endif
