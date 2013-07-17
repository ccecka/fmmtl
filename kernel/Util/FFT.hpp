#pragma once
//////////////////////
// FFTs and Helpers //
//////////////////////
// Link with -lfftw3 -lm

#include <fftw3.h>

#include <complex>
typedef std::complex<double> complex;

/** Scale the sequence elements by a constant */
template <class T, class Ts>
inline T* scale(T* a, int N, Ts scale) {
  for (int k = 0; k < N; ++k)
    a[k] *= scale;
  return a;
}

/** Truncates the sequence (stride S) from N0 elements to N elements
 * by removing the elements in the center of the sequence
 * @pre N0 > N */
template <class T>
inline T* f_trunc(T* a, int N0, int N, int S = 1) {
  int M = N/2;                       // Note integer division
  int k = M + (N & 1);
  int index = N0 - M;
  for ( ; k < N; ++k, ++index) {
    a[k*S] = a[index*S];
    a[index*S] = T(0);
  }
  // Zero out anything that was missed
  for ( ; k < N0 - M; ++k)
    a[k*S] = T(0);
  return a;
}

/** Extend the sequence (stride S) from N0 elements to N elements
 * by inserting zero elements into the center of the sequence
 * @pre N0 < N */
template <class T>
inline T* f_extend(T* a, int N0, int N, int S = 1) {
  int end = N0/2 - !(N0 & 1);       // Note integer division
  int index = N0 - 1;
  int q = N - 1;
  for ( ; index > end; --q, --index) {
    a[q*S] = a[index*S];
    a[index*S] = T(0);
  }
  // Zero out anything that was missed
  for ( ; q > N0 - 1; --q)
    a[q*S] = T(0);
  return a;
}

/** Truncate or extend the sequence (stride S) from N0 elements to N elements
 * by inserting or removing elements in the center of the range.
 * This method also scales the result with N/N0 */
template <class T, class Ts = double>
inline T* f_smooth(T* a, int N0, int N, int S = 1) {
  if (N0 > N) {
    return scale(f_trunc(a,N0,N,S), N, double(N)/N0);
  } else if (N0 < N) {
    return f_extend(scale(a,N0,double(N)/N0), N0, N, S);
  } // else they're equal and do nothing
  return a;
}

/** Truncate or extend the sequence (stride S) from N0 elements to N elements
 * by inserting or removing elements in the center of the range.
 * This method also scales the result with N/N0 */
template <class T>
inline T* f_cut(T* a, int N0, int N, int S = 1) {
  if (N0 > N) {
    return f_trunc(a, N0, N, S);
  } else if (N0 < N) {
    return f_extend(a, N0, N, S);
  } // else they're equal and do nothing
  return a;
}

// Perform a (slow) forward fft
inline complex* fft(complex* in, complex* out, int N)  {
  fftw_complex* a = reinterpret_cast<fftw_complex*>(in);
  fftw_complex* b = reinterpret_cast<fftw_complex*>(out);
  fftw_plan p = fftw_plan_dft_1d(N, a, b, FFTW_FORWARD, FFTW_ESTIMATE);
  fftw_execute(p);
  fftw_destroy_plan(p);
  return out;
}

// Perform a (slow) in-place forward fft
inline complex* fft(complex* in, int N) {
  return fft(in,in,N);
}

// Perform a (slow) backward fft
inline complex* ifft(complex* in, complex* out, int N) {
  fftw_complex* a = reinterpret_cast<fftw_complex*>(in);
  fftw_complex* b = reinterpret_cast<fftw_complex*>(out);
  fftw_plan p = fftw_plan_dft_1d( N, a, b, FFTW_BACKWARD, FFTW_ESTIMATE );
  fftw_execute(p);
  fftw_destroy_plan(p);
  scale(out, N, 1.0/N);
  return out;
}

// Perform a (slow) in-place backward fft
inline complex* ifft(complex* in, int N) {
  return ifft(in,in,N);
}

// Perform a (slow) in-place fft interpolation
inline complex* fftinterp(complex* in, int N0, int N)
{
  if (N0 == N)
    return in;

  fftw_complex* a = reinterpret_cast<fftw_complex*>(in);
  fftw_plan pFFT = fftw_plan_dft_1d(N0, a, a, FFTW_FORWARD, FFTW_ESTIMATE);
  fftw_plan pIFFT = fftw_plan_dft_1d(N, a, a, FFTW_BACKWARD, FFTW_ESTIMATE);

  fftw_execute(pFFT);
  if (N0 > N)       scale(f_trunc(in, N0, N), N, 1.0/N0);
  else if (N0 < N)  f_extend(scale(in, N0, 1.0/N0), N0, N);
  fftw_execute( pIFFT );

  fftw_destroy_plan( pFFT );
  fftw_destroy_plan( pIFFT );
  return in;
}
