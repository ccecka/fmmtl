#pragma once

#include <complex>
#include <cmath>

/** Spherical to cartesian coordinates */
template <typename real_type, typename point_type>
inline point_type sph2cart(real_type rho, real_type theta, real_type phi,
                           const point_type& s) {
  real_type st = std::sin(theta);
  real_type ct = std::cos(theta);
  real_type sp = std::sin(phi);
  real_type cp = std::cos(phi);
  return point_type(s[0]*st*cp + s[1]*ct*cp/rho - s[2]*sp/(rho*st),
                    s[0]*st*sp + s[1]*ct*sp/rho + s[2]*cp/(rho*st),
                    s[0]*ct    - s[1]*st/rho);
}

/** Cartesian to spherical coordinates */
template <typename real_type, typename point_type>
inline void cart2sph(real_type& r, real_type& theta, real_type& phi,
                     const point_type& x) {
  r = norm(x);
  theta = std::acos(x[2] / (r + 1e-100));
  phi = std::atan2(x[1], x[0]);
}



/** Computes
 * Y[n*(n+1)+m] = A_n^m rho^n Y_n^m(theta, phi)
 *              = (-1)^n/sqrt((n+m)!(n-m)!) rho^n Y_n^m(theta, phi)
 *              = (-1)^n/(n+m)! rho^n P_n^m(cos theta) exp(i m phi)
 * for all 0 <= n < P and all -n <= m <= n.
 *
 * Note that this uses the definition
 *
 * ??? How different from evalLocal computes?
 *
 * Note these are not the spherical harmonics, but are the spherical
 * harmonics with the prefactor (often denoted A_n^m) included. These are useful
 * for computing multipole and local expansions in an FMM.
 */
template <typename real>
void evalMultipole(real rho, real theta, real phi, int P,
                   std::complex<real>* Y, std::complex<real>* dY = nullptr) {
  typedef std::complex<real> complex;
  using std::cos;
  using std::sin;
  const real    ct = cos(theta);
  const real    st = sin(theta);
  const complex ei = complex(cos(phi), sin(phi)); // exp(i * phi)
  real    Pmm = 1;                                // Init Legendre P00(ct)
  real   rhom = 1;                                // Init (-1)^n rho^n / (n+m)!
  complex eim = 1;                                // Init exp(i*m*phi)
  int m = 0;
  while (true) {
    // n == m
    int npn = m * m + 2 * m;                      //  Index of Ynm for m > 0
    int nmn = m * m;                              //  Index of Ynm for m < 0
    Y[npn] = rhom * Pmm * eim;                    //  Ynm for m > 0
    Y[nmn] = std::conj(Y[npn]);                   //  Conj for m < 0
    if (dY)
      dY[npn] = m*ct/st * Y[npn];                 // theta derivative

    real Pn1m = Pmm;                              //  P_m^m
    real Pnm  = ct * (2*m+1) * Pmm;               //  P_{m+1}^{m}(x) = x (2m+1) Pmm
    real rhon = rhom;                             //  (-1)^m rho^m / (2m)!
    for (int n = m+1; n < P; ++n) {
      rhon *= -rho / (n + m);                     //   (-1)^n rho^n / (n+m)!
      int npm = n * n + n + m;                    //   Index of Ynm for m > 0
      int nmm = n * n + n - m;                    //   Index of Ynm for m < 0

      Y[npm] = rhon * Pnm * eim;                  //   Ynm for m > 0
      Y[nmm] = std::conj(Y[npm]);                 //   Conj for m < 0
      if (dY)
        dY[npm] = (n*ct - (n+m)*Pn1m/Pnm)/st * Y[npm];  // theta derivative

      real Pn2m = Pn1m;                           //   P_{n-1}^m
      Pn1m = Pnm;                                 //   P_n^m
      Pnm = (ct*(2*n+1)*Pn1m-(n+m)*Pn2m)/(n-m+1); //   P_{n+1}^m recurrence
    }                                             //  End loop over n in Ynm

    ++m;                                          // Increment m
    if (m == P) return;

    rhom *= -rho / (2*m*(2*m-1));                 //  (-1)^m rho^m / (2m)!
    Pmm = -st * (2*m-1) * Pmm;                    //  P_{m+1}^{m+1} recurrence
    eim *= ei;                                    //  exp(i*m*phi)
  }                                               // End loop over m in Ynm
}


/** Computes
 * Y[n*(n+1)+m] = A_n^m i^-|m| (-rho)^n Y_n^m(theta, phi)
 *              = (-1)^n/sqrt((n+m)!(n-m)!) i^-|m| (-rho)^n Y_n^m(theta, phi)
 *              = rho^n/(n+m)! P_n^m(cos theta) exp(i m phi) i^-|m|
 * for all 0 <= n < P and all -n <= m <= n.
 *
 * Note that this uses the definition...
 *
 * Note these are not the spherical harmonics, but are the spherical
 * harmonics with the prefactor (often denoted A_n^m) included. These are useful
 * for computing multipole and local expansions in an FMM.
 */
template <typename real>
void evalZ(real rho, real theta, real phi, int P,
           std::complex<real>* Y, std::complex<real>* dY = nullptr) {
  typedef std::complex<real> complex;
  using std::cos;
  using std::sin;
  for (int i = 0; i < P*P; ++i) Y[i] = std::numeric_limits<double>::quiet_NaN();

  const real    ct = cos(theta);
  const real    st = sin(theta);
  const complex ei = complex(sin(phi),-cos(phi)); // exp(i*phi) i^-1

  real    Pmm = 1;                                // Init Legendre P00(ct)
  real   rhom = 1;                                // Init rho^n / (n+m)!
  complex eim = 1;                                // Init exp(i*m*phi) i^-m
  int m = 0;
  while (true) {
    // n == m
    int npn = m*(m+1) + m;                        //  Index of Ynm for m > 0
    Y[npn] = rhom * Pmm * eim;                    //  Ynm for m > 0
    if (dY)
      dY[npn] = m*ct/st * Y[npn];                 // theta derivative

    // n == m+1
    int n = m + 1;
    if (n == P) return;                           // Done! m == P-1

    real Pnm  = ct * (2*m+1) * Pmm;               //  P_{m+1}^m(x) = x(2m+1)Pmm
    real rhon = rhom * rho / (n+m);               //  rho^n / (n+m)!
    int npm = n*(n+1) + m;                        //   Index of Ynm for m > 0
    Y[npm] = rhon * Pnm * eim;                    // Ynm for m > 0
    if (dY)
      dY[npm] = (n*ct - (n+m)*Pmm/Pnm)/st * Y[npm];  // theta derivative

    // m+1 < n < P
    real Pn1m = Pmm;                              //  P_{n-1}^m
    while (++n != P) {
      real Pn2m = Pn1m;                           //   P_{n-2}^m
      Pn1m = Pnm;                                 //   P_{n-1}^m
      Pnm = (ct*(2*n-1)*Pn1m-(n+m-1)*Pn2m)/(n-m); //   P_n^m recurrence
      rhon *= rho / (n + m);                      //   rho^n / (n+m)!

      int npm = n*(n+1) + m;                      //   Index of Ynm for m > 0
      Y[npm] = rhon * Pnm * eim;                  //   Ynm for m > 0
      if (dY)
        dY[npm] = (n*ct - (n+m)*Pn1m/Pnm)/st * Y[npm];  // theta derivative
    }

    ++m;                                          // Increment m

    rhom *= rho / (2*m*(2*m-1));                  //  rho^m / (2m)!
    Pmm *= -st * (2*m-1);                         //  P_{m+1}^{m+1} recurrence
    eim *= ei;                                    //  exp(i*m*phi) i^-m
  }                                               // End loop over m in Ynm
}






//! Evaluate singular harmonics \f$ r^{-n-1} Y_n^m \f$
template <typename real>
void evalLocal(real rho, real theta, real phi, int P,
               std::complex<real>* Ynm) {
  typedef std::complex<real> complex;
  const real    ct = cos(theta);
  const real    st = sin(theta);
  const complex ei = complex(cos(phi), sin(phi)); // exp(i * phi)
  real Pmm = 1;                                                // Initialize Legendre polynomial Pmm
  real invR = -1.0 / rho;                                     // - 1 / rho
  real rhom = -invR;                                          // Initialize rho^(-m-1)
  complex eim = 1.0;                                          // Initialize exp(i * m * phi)
  for (int m=0; m<P; m++) {                                     // Loop over m in Ynm
    int npn = m * m + 2 * m;                                    //  Index of Ynm for m > 0
    int nmn = m * m;                                            //  Index of Ynm for m < 0
    Ynm[npn] = rhom * Pmm * eim;                                  //  rho^(-m-1) * Ynm for m > 0
    Ynm[nmn] = std::conj(Ynm[npn]);                             //  Use conjugate relation for m < 0

    real p1 = Pmm;                                              //  Pnm-1
    real p = ct * (2 * m + 1) * p1;                                   //  Pnm using recurrence relation
    rhom *= invR;                                               //  rho^(-m-1)
    real rhon = rhom;                                         //  rho^(-n-1)
    for (int n=m+1; n<P; n++) {                                 //  Loop over n in Ynm
      int npm = n * n + n + m;                                  //   Index of Ynm for m > 0
      int nmm = n * n + n - m;                                  //   Index of Ynm for m < 0
      Ynm[npm] = rhon * p * eim;                                //   rho^n * Ynm for m > 0
      Ynm[nmm] = std::conj(Ynm[npm]);                           //   Use conjugate relation for m < 0
      real p2 = p1;                                           //   Pnm-2
      p1 = p;                                                   //   Pnm-1
      p = (ct * (2 * n + 1) * p1 - (n + m) * p2) / (n - m + 1);  //   Pnm using recurrence relation
      rhon *= invR * (n - m + 1);                               //   rho^(-n-1)
    }
                                                        //  End loop over n in Ynm
    Pmm = -Pmm * (2*m+1) * st;                                        //  Pn
    eim *= ei;                                                  //  Update exp(i * m * phi)
  }                                                             // End loop over m in Ynm
}


/** Computes the function
 * Y[n*(n+1)+m] = i^|m| / A_n^m rho^{-n-1} Y_n^m(theta, phi)
 *              = i^|m| (-1)^n sqrt((n+m)!(n-m)!) rho^{-n-1} Y_n^m(theta, phi)
 *              = (-1)^n rho^{-n-1} (n-|m|)! P_n^|m|(cos theta) exp(i m phi) i^|m|
 * for all 0 <= n < P and all -n <= m <= n.
 *
 * Note that this uses the definition...
 *
 * Note these are not the spherical harmonics, but are the spherical
 * harmonics with the prefactor (often denoted A_n^m) included. These are useful
 * for computing multipole and local expansions in an FMM.
 */
template <typename real>
void evalW(real rho, real theta, real phi, int P,
           std::complex<real>* Y, std::complex<real>* dY = nullptr) {
  typedef std::complex<real> complex;
  using std::cos;
  using std::sin;
  for (int i = 0; i < P*P; ++i) Y[i] = std::numeric_limits<double>::quiet_NaN();

  rho = 1 / rho;

  const real    ct = cos(theta);
  const real    st = sin(theta);
  const complex ei = complex(-sin(phi),cos(phi)); // exp(i*phi) i

  real    Pmm = 1;                                // Init Legendre P00(ct)
  real   rhom = rho;                              // Init (-1)^n rho^{-n-1} (n-m)!
  complex eim = 1;                                // Init exp(i*m*phi) i^-m
  int m = 0;
  while (true) {
    // n == m
    int npn = m*(m+1) + m;                        //  Index of Ynm for m > 0
    Y[npn] = rhom * Pmm * eim;                    //  Ynm for m > 0
    if (dY)
      dY[npn] = m*ct/st * Y[npn];                 // theta derivative

    // n == m+1
    int n = m+1;
    if (n == P) return;                           // Done! m == P-1

    real Pnm  = ct * (2*m+1) * Pmm;               //  P_{m+1}^m(x) = x(2m+1)Pmm
    real rhon = rhom * -rho;                      //  (-1)^n rho^{-n-1} (n-m)!
    int npm = n*(n+1) + m;                        //   Index of Ynm for m > 0
    Y[npm] = rhon * Pnm * eim;                    // Ynm for m > 0
    if (dY)
      dY[npm] = (n*ct - (n+m)*Pmm/Pnm)/st * Y[npm];  // theta derivative

    // m+1 < n < P
    real Pn1m = Pmm;                              //  P_{n-1}^m
    while (++n != P) {
      real Pn2m = Pn1m;                           //   P_{n-2}^m
      Pn1m = Pnm;                                 //   P_{n-1}^m
      Pnm = (ct*(2*n-1)*Pn1m-(n+m-1)*Pn2m)/(n-m); //   P_n^m recurrence
      rhon *= -rho * (n - m);                     //   (-1)^n rho^{-n-1} (n-m)!

      int npm = n*(n+1) + m;                      //   Index of Ynm for m > 0
      Y[npm] = rhon * Pnm * eim;                  //   Ynm for m > 0
      if (dY)
        dY[npm] = (n*ct - (n+m)*Pn1m/Pnm)/st * Y[npm];  // theta derivative
    }

    ++m;                                          // Increment m

    rhom *= -rho;                                 //  (-1)^m rho^{-m-1} (n-m)!
    Pmm *= -st * (2*m-1);                         //  P_{m+1}^{m+1} recurrence
    eim *= ei;                                    //  exp(i*m*phi) i^m
  }                                               // End loop over m in Ynm
}
