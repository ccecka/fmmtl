#include <complex>
#include <cmath>

#if 0
//! Evaluate solid harmonics \f$ a r^n Y_n^m \f$
/** Computes Ynm[n*(n+1)+m] = a rho^n Y_n^m(theta, phi)
 * with
 * Y_n^m(theta,phi) = (-1)^m sqrt((n-|m|)!/(n+|m|)!) P_n^|m|(cos theta) exp(i m phi)
 * PLUS THE PREFACTOR AND AND STUFF
 */
template <typename T>
inline void
generate_Y_2(T rho, T theta, T phi, T a, unsigned P,
             std::complex<T>* Ynm, std::complex<T>* YnmTheta) {
  typedef T real;
  typedef std::complex<T> complex;
  const real ct = std::cos(theta);                        // ct = cos(theta)
  const real st = std::sin(theta);                        // st = sin(theta)
  const complex ei = std::exp(complex(0,phi));            // exp(i*phi)
  real pn = 1;                                            // Init Legendre Pn
  real rhom = a;                                          // Init rho^m
  complex eim = 1;                                        // Init exp(i*m*phi)
  for (unsigned m = 0; m != P; ++m) {                     // Loop m in Ynm
    int npn = m * m + 2 * m;                              //  Index for m > 0
    int nmn = m * m;                                      //  Index for m < 0
    Ynm[npn] = (rhom * pn) * eim;                         //  rho^m * Ynm
    Ynm[nmn] = std::conj(Ynm[npn]);                       //  Conj for m < 0
    YnmTheta[npn] = (m*ct/st) * Ynm[npn];                 //  Theta derivative
    YnmTheta[nmn] = std::conj(YnmTheta[npn]);             //  Conj for m < 0

    real p1 = pn;                                         //  Pnm-1
    real p = ct * (2*m+1) * pn;                           //  Pnm recurrence
    real rhon = rhom;                                     //  rho^n
    for (unsigned n = m+1; n != P; ++n) {                 //  Loop n in Ynm
      rhon *= -rho / (n+m);                               //   Update factorial
      int npm = n * n + n + m;                            //   Index for m > 0
      int nmm = n * n + n - m;                            //   Index for m < 0
      Ynm[npm] = (rhon * p) * eim;                        //   rho^n * Ynm
      Ynm[nmm] = std::conj(Ynm[npm]);                     //   Conj for m < 0
      YnmTheta[npm] = ((n*ct - (n+m)*p1/p)/st) * Ynm[npm];//   Theta derivative
      YnmTheta[nmm] = std::conj(YnmTheta[npm]);           //   Conj for m < 0

      real p2 = p1;                                       //   Pnm-2
      p1 = p;                                             //   Pnm-1
      p = ((2*n+1)*ct*p1 - (n+m)*p2) / (n-m+1);           //   Pnm recurrence
    }                                                     //  End loop n in Ynm
    rhom *= -rho / ((2*m+2)*(2*m+1));                     //  Update factorial
    pn *= -(2*m+1) * st;                                  //  Pn recurrence
    eim *= ei;                                            //  exp(i*m*phi)
  }                                                       // End loop m in Ynm
}
#endif

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





//! Evaluate solid harmonics \f$ r^n Y_{n}^{m} \f$
/** Computes Ynm[n*(n+1)+m] = rho^n Y_n^m(theta, phi)
 * with
 * Y_n^m(theta,phi) = (-1)^m sqrt((n-|m|)!/(n+|m|)!) P_n^|m|(cos theta) exp(i m phi)
 * PLUS THE PREFACTOR AND AND STUFF
 */
template <typename T>
void evalMultipole(T rho, T theta, T phi, int P,
                   std::complex<T>* Ynm) {
  typedef T real_type;
  typedef std::complex<T> complex_type;
  real_type x = std::cos(theta);                                   // x = cos(theta)
  real_type y = std::sin(theta);                                   // y = sin(theta)
  real_type pn = 1;                                                // Initialize Legendre polynomial Pn
  real_type rhom = 1;                                              // Initialize rho^m
  complex_type ei = std::exp(complex_type(0,phi));                   // exp(i * phi)
  complex_type eim = 1.0;                                          // Initialize exp(i * m * phi)
  for (int m=0; m<P; m++) {                                     // Loop over m in Ynm
    real_type p = pn;                                              //  Associated Legendre polynomial Pnm
    int npn = m * m + 2 * m;                                    //  Index of Ynm for m > 0
    int nmn = m * m;                                            //  Index of Ynm for m < 0
    Ynm[npn] = rhom * p * eim;                                  //  rho^m * Ynm for m > 0
    Ynm[nmn] = std::conj(Ynm[npn]);                             //  Use conjugate relation for m < 0
    real_type p1 = p;                                              //  Pnm-1
    p = x * (2 * m + 1) * p1;                                   //  Pnm using recurrence relation
    rhom *= rho;                                                //  rho^m
    real_type rhon = rhom;                                         //  rho^n
    for (int n=m+1; n<P; n++) {                                 //  Loop over n in Ynm
      int npm = n * n + n + m;                                  //   Index of Ynm for m > 0
      int nmm = n * n + n - m;                                  //   Index of Ynm for m < 0
      rhon /= -(n + m);                                         //   Update factorial
      Ynm[npm] = rhon * p * eim;                                //   rho^n * Ynm
      Ynm[nmm] = std::conj(Ynm[npm]);                           //   Use conjugate relation for m < 0
      real_type p2 = p1;                                           //   Pnm-2
      p1 = p;                                                   //   Pnm-1
      p = (x * (2 * n + 1) * p1 - (n + m) * p2) / (n - m + 1);  //   Pnm using recurrence relation
      rhon *= rho;                                              //   Update rho^n
    }                                                           //  End loop over n in Ynm
    rhom /= -(2 * m + 2) * (2 * m + 1);                         //  Update factorial
    pn = -pn * (2*m+1) * y;                                        //  Pn
    eim *= ei;                                                  //  Update exp(i * m * phi)
  }                                                             // End loop over m in Ynm
}



//! Evaluate solid harmonics \f$ r^n Y_{n}^{m} \f$
template <typename T>
void evalMultipole(T rho, T theta, T phi, int P,
                   std::complex<T>* Ynm, std::complex<T>* YnmTheta) {
  typedef T real_type;
  typedef std::complex<T> complex_type;
  real_type x = std::cos(theta);                                   // x = cos(theta)
  real_type y = std::sin(theta);                                   // y = sin(theta)
  real_type pn = 1;                                                // Initialize Legendre polynomial Pn
  real_type rhom = 1;                                              // Initialize rho^m
  complex_type ei = std::exp(complex_type(0,phi));                   // exp(i * phi)
  complex_type eim = 1.0;                                          // Initialize exp(i * m * phi)
  for (int m=0; m<P; m++) {                                     // Loop over m in Ynm
    real_type p = pn;                                              //  Associated Legendre polynomial Pnm
    int npn = m * m + 2 * m;                                    //  Index of Ynm for m > 0
    int nmn = m * m;                                            //  Index of Ynm for m < 0
    Ynm[npn] = rhom * p * eim;                                  //  rho^m * Ynm for m > 0
    Ynm[nmn] = std::conj(Ynm[npn]);                             //  Use conjugate relation for m < 0
    real_type p1 = p;                                              //  Pnm-1
    p = x * (2 * m + 1) * p1;                                   //  Pnm using recurrence relation
    YnmTheta[npn] = rhom * (p - (m + 1) * x * p1) / y * eim;    // theta derivative of r^n * Ynm
    rhom *= rho;                                                //  rho^m
    real_type rhon = rhom;                                         //  rho^n
    for (int n=m+1; n<P; n++) {                                 //  Loop over n in Ynm
      int npm = n * n + n + m;                                  //   Index of Ynm for m > 0
      int nmm = n * n + n - m;                                  //   Index of Ynm for m < 0
      rhon /= -(n + m);                                         //   Update factorial
      Ynm[npm] = rhon * p * eim;                                //   rho^n * Ynm
      Ynm[nmm] = std::conj(Ynm[npm]);                           //   Use conjugate relation for m < 0
      real_type p2 = p1;                                           //   Pnm-2
      p1 = p;                                                   //   Pnm-1
      p = (x * (2 * n + 1) * p1 - (n + m) * p2) / (n - m + 1);  //   Pnm using recurrence relation
      YnmTheta[npm] = rhon * ((n - m + 1) * p - (n + 1) * x * p1) / y * eim;// theta derivative
      rhon *= rho;                                              //   Update rho^n
    }                                                           //  End loop over n in Ynm
    rhom /= -(2 * m + 2) * (2 * m + 1);                         //  Update factorial
    pn = -pn * (2*m+1) * y;                                        //  Pn
    eim *= ei;                                                  //  Update exp(i * m * phi)
  }                                                             // End loop over m in Ynm
}

//! Evaluate singular harmonics \f$ r^{-n-1} Y_n^m \f$
template <typename T>
void evalLocal(T rho, T theta, T phi, int P,
               std::complex<T>* Ynm) {
  typedef T real_type;
  typedef std::complex<T> complex_type;
  real_type x = std::cos(theta);                                   // x = cos(theta)
  real_type y = std::sin(theta);                                   // y = sin(theta)
  real_type fact = 1;                                              // Initialize 2 * m + 1
  real_type pn = 1;                                                // Initialize Legendre polynomial Pn
  real_type invR = -1.0 / rho;                                     // - 1 / rho
  real_type rhom = -invR;                                          // Initialize rho^(-m-1)
  complex_type ei = std::exp(complex_type(0,phi));                            // exp(i * phi)
  complex_type eim = 1.0;                                          // Initialize exp(i * m * phi)
  for (int m=0; m<P; m++) {                                     // Loop over m in Ynm
    real_type p = pn;                                              //  Associated Legendre polynomial Pnm
    int npn = m * m + 2 * m;                                    //  Index of Ynm for m > 0
    int nmn = m * m;                                            //  Index of Ynm for m < 0
    Ynm[npn] = rhom * p * eim;                                  //  rho^(-m-1) * Ynm for m > 0
    Ynm[nmn] = std::conj(Ynm[npn]);                             //  Use conjugate relation for m < 0
    real_type p1 = p;                                              //  Pnm-1
    p = x * (2 * m + 1) * p1;                                   //  Pnm using recurrence relation
    rhom *= invR;                                               //  rho^(-m-1)
    real_type rhon = rhom;                                         //  rho^(-n-1)
    for (int n=m+1; n<P; n++) {                                 //  Loop over n in Ynm
      int npm = n * n + n + m;                                  //   Index of Ynm for m > 0
      int nmm = n * n + n - m;                                  //   Index of Ynm for m < 0
      Ynm[npm] = rhon * p * eim;                                //   rho^n * Ynm for m > 0
      Ynm[nmm] = std::conj(Ynm[npm]);                           //   Use conjugate relation for m < 0
      real_type p2 = p1;                                           //   Pnm-2
      p1 = p;                                                   //   Pnm-1
      p = (x * (2 * n + 1) * p1 - (n + m) * p2) / (n - m + 1);  //   Pnm using recurrence relation
      rhon *= invR * (n - m + 1);                               //   rho^(-n-1)
    }                                                           //  End loop over n in Ynm
    pn = -pn * fact * y;                                        //  Pn
    fact += 2;                                                  //  2 * m + 1
    eim *= ei;                                                  //  Update exp(i * m * phi)
  }                                                             // End loop over m in Ynm
}
