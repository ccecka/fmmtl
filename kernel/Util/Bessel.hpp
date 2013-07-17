#pragma once
//////////////////
// GSL Wrappers //
//////////////////
// Link with -lgsl -lm

#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_sf_legendre.h>

#include <gsl/gsl_errno.h>

#include <complex>

// Override GSL underflow error exit
struct _GSLError_ {
  _GSLError_() { gsl_set_error_handler(&_GSLError_::Handler); }
  static void Handler(const char* msg, const char* file, int line, int err) {
    printf("GSLError %d in %s at %d : %s\n",err,file,line,msg);
    if (err != GSL_EUNDRFLW)
      exit(1);
  }
};
// Re-define GSL default error handler when loading the library
static _GSLError_ __GSLError__;


/** Cylindrical Bessel function J_n(x)
 * Returns 0 in the case of underflow */
inline double bessel_J(int n, double x) {
  gsl_sf_result result;
  int status = gsl_sf_bessel_Jn_e(n,x,&result);
  if (status == GSL_EUNDRFLW) return 0;
  return result.val;
}

/** Spherical Bessel function j_n(x)
 * Returns 0 in the case of underflow */
inline double bessel_j(int n, double x) {
  gsl_sf_result result;
  int status = gsl_sf_bessel_jl_e(n,x,&result);
  if (status == GSL_EUNDRFLW) return 0;
  return result.val;
}

/** Cylindral Bessel function Y_n(x)
 * Returns 0 in the case of underflow */
inline double bessel_Y(int n, double x) {
  gsl_sf_result result;
  int status = gsl_sf_bessel_Yn_e(n,x,&result);
  if (status == GSL_EUNDRFLW) return 0;
  return result.val;
}

/** Spherical Bessel function y_n(x)
 * Returns 0 in the case of underflow */
inline double bessel_y(int n, double x) {
  gsl_sf_result result;
  int status = gsl_sf_bessel_yl_e(n,x,&result);
  if (status == GSL_EUNDRFLW) return 0;
  return result.val;
}

/** Spherical Hankel function h_n(x) */
inline std::complex<double> bessel_h(int n, double x) {
  return std::complex<double>(bessel_j(n,x), bessel_y(n,x));
}

/** Legendre polynomial P_n(x)
 * Returns 0 in the case of underflow */
inline double legendre_P(int n, double x) {
  gsl_sf_result result;
  int status = gsl_sf_legendre_Pl_e(n,x,&result);
  if (status == GSL_EUNDRFLW) return 0;
  return result.val;
}
