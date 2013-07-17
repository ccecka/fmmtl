#pragma once

#include "FFT.hpp"

#include <vector>
#include <complex>
#include <algorithm>

// A reterpolator between two quadratures using FFT interpolation
// Optimized to reduce the number of FFTs
// Generalized to handle any quadrature
class FFTInterp_FFT2 {
  S2Quadrature* q1;
  S2Quadrature* q2;

  // q1 Data
  unsigned N_rows1, N_theta1, N_phi1;

  // q2 Data
  unsigned N_rows2, N_theta2, N_phi2;

  // Intermediate Data
  unsigned N_rows, N_phi, N_phi0;

  // Intermediate Storage, size = N_rows x N_phi
  std::vector<complex> T;
  // Intermediate Storage, size = max(N_theta1, N_theta2)
  std::vector<complex> t;

  // First stage phi (row) FFTs
  std::vector<fftw_plan> rFFT1;

  // Second stage theta (column) FFTs
  fftw_plan tFFT;
  fftw_plan tIFFT;

  // Third stage phi (row) IFFTs
  std::vector<fftw_plan> rIFFT2;

 public:

  // Constructor
  FFTInterp_FFT2(S2Quadrature* q1_, S2Quadrature* q2_)
      : q1(q1_), q2(q2_),
        N_rows1(q1->n_row()), N_theta1(q1->n_theta()), N_phi1(q1->n_phi()),
        N_rows2(q2->n_row()), N_theta2(q2->n_theta()), N_phi2(q2->n_phi()),
        N_rows(std::max(N_rows1,N_rows2)),
        N_phi(std::max(N_phi1,N_phi2)),
        N_phi0(std::min(N_phi1,N_phi2)),
        T(N_phi * N_rows),
        t(std::max(N_theta1,N_theta2))
  {
    // Define an phi1 FFT for each row
    rFFT1 = std::vector<fftw_plan>(N_rows1);
    for (unsigned n = 0; n < N_rows1; ++n) {
      fftw_complex* Tn = reinterpret_cast<fftw_complex*>(T.data()+n*N_phi);
      unsigned Nphi_n = q1->n_phi(n);
      rFFT1[n] = fftw_plan_dft_1d(Nphi_n,Tn,Tn,FFTW_FORWARD,FFTW_MEASURE);
    }

    // Define a theta FFT and IFFT (could also do these as blocks...)
    fftw_complex* t0 = reinterpret_cast<fftw_complex*>(t.data());
    tFFT  = fftw_plan_dft_1d(N_theta1, t0, t0, FFTW_FORWARD,  FFTW_MEASURE);
    tIFFT = fftw_plan_dft_1d(N_theta2, t0, t0, FFTW_BACKWARD, FFTW_MEASURE);

    // Define an phi2 IFFT for each row
    rIFFT2 = std::vector<fftw_plan>(N_rows2);
    for (unsigned n = 0; n < N_rows2; ++n) {
      fftw_complex* Tn = reinterpret_cast<fftw_complex*>(T.data()+n*N_phi);
      unsigned Nphi_n = q2->n_phi(n);
      rIFFT2[n] = fftw_plan_dft_1d(Nphi_n,Tn,Tn,FFTW_BACKWARD,FFTW_MEASURE);
    }
  }

  // Destructor
  ~FFTInterp_FFT2() {
    std::for_each(rFFT1.begin(), rFFT1.end(), fftw_destroy_plan);
    fftw_destroy_plan(tFFT);
    fftw_destroy_plan(tIFFT);
    std::for_each(rIFFT2.begin(), rIFFT2.end(), fftw_destroy_plan);
  }

  inline void apply(const S2Function& A, S2Function& B)
  {
    assert(A.quad == q1 && B.quad == q2);

    // Copy from A and interpolate each row
    for (unsigned n = 0; n < N_rows1; ++n) {
      unsigned N_phi_n = q1->n_phi(n);
      complex* Tn = T.data() + n*N_phi;

      // Copy into the nth row of T
      std::copy(A.row_begin(n), A.row_end(n), Tn);
      // FFT the nth row of T (FFT in phi)
      fftw_execute(rFFT1[n]);
      // Add/remove the high frequencies and scale for the upcoming interpolations
      scale(Tn, N_phi_n, 1.0/(N_phi_n*N_theta1));
      f_cut(Tn, N_phi_n, N_phi0);
    }

    // Interpolate from N_theta1 to N_theta2
    for (unsigned m = 0; m < N_phi0; ++m) {
      double one = (m & 1) ? -1.0 : 1.0;

      // Unwrap into t
      t[0] = T[m + 0*N_phi];
      for (unsigned n = 1; n < N_rows1; ++n)
        t[N_theta1-n] = one * (t[n] = T[m + n*N_phi]);

      // Interpolate from N_theta1 to N_theta2
      fftw_execute(tFFT);                          // FFT_theta
      f_cut(t.data(), N_theta1, N_theta2);         // Add/remove high frequencies
      fftw_execute(tIFFT);                         // IFFT_theta

      // Unwrap back into T
      for (unsigned n = 0; n < N_rows2; ++n)
        T[m + n*N_phi] = t[n];
    }

    // Interpolate each row and copy into B
    for (unsigned n = 0; n < N_rows2; ++n) {
      unsigned N_phi_n = q2->n_phi(n);
      complex* Tn = T.data() + n*N_phi;

      // Add/remove the high frequencies
      f_cut(Tn, N_phi0, N_phi_n);
      // IFFT the nth row of T (IFFT in phi)
      fftw_execute(rIFFT2[n]);

      // Copy nth row of T into B
      std::copy(Tn, Tn + N_phi_n, B.row_begin(n));
    }
  }

 private:
  // Disable Copy and Assignment
  FFTInterp_FFT2(const FFTInterp_FFT2&) {}
  void operator=(const FFTInterp_FFT2&) {}
};

typedef FFTInterp_FFT2 Reterpolator;
