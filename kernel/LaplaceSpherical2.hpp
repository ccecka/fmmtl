#pragma once
/** @file LaplaceSpherical.hpp
 * @brief Implements the Laplace kernel with spherical expansions.
 *
 * K(t,s) = 1 / |s-t|        // Laplace potential
 * K(t,s) = (s-t) / |s-t|^3  // Laplace force
 */

#include <complex>
#include <cmath>
#include <cassert>

#include "Laplace.kern"

#include "fmmtl/Expansion.hpp"
// Use a library-defined Vector class that supports multiple architectures
#include "fmmtl/numeric/Vec.hpp"

#define ODDEVEN(n) ((((n) & 1) == 1) ? -1 : 1)
#define IPOW2N(n) ((n >= 0) ? 1 : ODDEVEN(n))

class LaplaceSpherical2
    : public fmmtl::Expansion<LaplaceKernel, LaplaceSpherical2> {
 protected:
  typedef double real;
  typedef std::complex<real> complex;

  //! Expansion order
  int P;

 public:
  //! Point type
  typedef Vec<3,real> point_type;

  //! Multipole expansion type
  typedef std::vector<complex> multipole_type;
  //! Local expansion type
  typedef std::vector<complex> local_type;

  //! Default constructor -- use delegating constructor
  LaplaceSpherical2() : LaplaceSpherical2(5) {}

  //! Constructor
  LaplaceSpherical2(int _P)
      : P(_P) {
  }

  /** Initialize a multipole expansion with the size of a box at this level */
  void init_multipole(multipole_type& M,
                      const point_type&, unsigned) const {
    M = std::vector<complex>(P*(P+1)/2);
  }
  /** Initialize a local expansion with the size of a box at this level */
  void init_local(local_type& L,
                  const point_type&, unsigned) const {
    L = std::vector<complex>(P*(P+1)/2);
  }

  /** Kernel P2M operation
   * M += Op(s) * c where M is the multipole and s is the source
   *
   * @param[in] source The point source
   * @param[in] charge The source's corresponding charge
   * @param[in] center The center of the box containing the multipole expansion
   * @param[in,out] M The multipole expansion to accumulate into
   */
  void P2M(const source_type& source, const charge_type& charge,
           const point_type& center, multipole_type& M) const {
    complex Ynm[P*P], YnmTheta[P*P];
    point_type r = source - center;
    real rho, alpha, beta;
    cart2sph(rho, alpha, beta, r);
    evalMultipole(rho, alpha, beta, Ynm, YnmTheta);
    for (int n = 0; n < P; ++n) {
      for (int m = 0; m <= n; ++m) {
        int nm  = n*(n+1)   - m;
        int nms = n*(n+1)/2 + m;
        M[nms] += charge * Ynm[nm];
      }
    }
  }

  /** Kernel M2M operator
   * M_t += Op(M_s) where M_t is the target and M_s is the source
   *
   * @param[in] source The multipole source at the child level
   * @param[in,out] target The multipole target to accumulate into
   * @param[in] translation The vector from source to target
   * @pre Msource includes the influence of all points within its box
   */
  void M2M(const multipole_type& Msource,
           multipole_type& Mtarget,
           const point_type& translation) const {
    complex Ynm[P*P], YnmTheta[P*P];
    real rho, alpha, beta;
    cart2sph(rho, alpha, beta, translation);
    evalMultipole(rho, alpha, beta, Ynm, YnmTheta);
    for (int j = 0; j < P; ++j) {
      for (int k = 0; k <= j; ++k) {
        int jks = j*(j+1)/2 + k;
        complex M = 0;
        for (int n = 0; n <= j; ++n) {
          for (int m = std::max(-n,-j+k+n); m <= std::min(k-1,n); ++m) {
            int jnkms = (j-n)*(j-n+1)/2 + k - m;
            int nm    = n*(n+1) - m;
            M += Msource[jnkms] * Ynm[nm] * real(IPOW2N(m) * ODDEVEN(n));
          }
          for (int m=k; m<=std::min(n,j+k-n); m++) {
            int jnkms = (j - n) * (j - n + 1) / 2 - k + m;
            int nm    = n * n + n - m;
            M += std::conj(Msource[jnkms]) * Ynm[nm] * real(ODDEVEN(k+n+m));
          }
        }
        Mtarget[jks] += M;
      }
    }
  }

  /** Kernel M2L operation
   * L += Op(M)
   *
   * @param[in] Msource The multpole expansion source
   * @param[in,out] Ltarget The local expansion target
   * @param[in] translation The vector from source to target
   * @pre translation obeys the multipole-acceptance criteria
   * @pre Msource includes the influence of all points within its box
   */
  void M2L(const multipole_type& Msource,
           local_type& Ltarget,
           const point_type& translation) const {
    complex Ynmi[P*P];
    real rho, alpha, beta;
    cart2sph(rho, alpha, beta, translation);
    evalLocal(rho, alpha, beta, Ynmi);
    for (int j = 0; j < P; ++j) {
      real Cnm = ODDEVEN(j);
      for (int k = 0; k <= j; k++) {
        int jks = j*(j+1)/2 + k;
        complex L = 0;
        for (int n=0; n<P-j; ++n) {
          for (int m=-n; m<0; ++m) {
            int nms  = n*(n+1)/2 - m;
            int jnkm = (j+n)*(j+n+1) + m - k;
            L += std::conj(Msource[nms]) * Cnm * Ynmi[jnkm];
          }
          for (int m = 0; m <= n; ++m) {
            int nms  = n*(n+1)/2 + m;
            int jnkm = (j+n)*(j+n+1) + m - k;
            real Cnm2 = Cnm * ODDEVEN((k-m)*(k<m)+m);
            L += Msource[nms] * Cnm2 * Ynmi[jnkm];
          }
        }
        Ltarget[jks] += L;
      }
    }
  }

  /** Kernel L2L operator
   * L_t += Op(L_s) where L_t is the target and L_s is the source
   *
   * @param[in] source The local source at the parent level
   * @param[in,out] target The local target to accumulate into
   * @param[in] translation The vector from source to target
   * @pre Lsource includes the influence of all points outside its box
   */
  void L2L(const local_type& Lsource,
           local_type& Ltarget,
           const point_type& translation) const {
    complex Ynm[P*P], YnmTheta[P*P];
    real rho, alpha, beta;
    cart2sph(rho, alpha, beta, translation);
    evalMultipole(rho, alpha, beta, Ynm, YnmTheta);
    for (int j = 0; j < P; ++j) {
      for (int k = 0; k <= j; ++k) {
        int jks = j*(j+1)/2 + k;
        complex L = 0;
        for (int n = j; n < P; ++n) {
          for (int m = j+k-n; m < 0; ++m) {
            int jnkm = (n-j)*(n-j+1) + m - k;
            int nms  = n*(n+1)/2 - m;
            L += std::conj(Lsource[nms]) * Ynm[jnkm] * real(ODDEVEN(k));
          }
          for (int m = 0; m <= n; ++m) {
            if (n-j >= abs(m-k)) {
              int jnkm = (n-j)*(n-j+1) + m - k;
              int nms  = n*(n+1)/2 + m;
              L += Lsource[nms] * Ynm[jnkm] * real(ODDEVEN((m-k)*(m<k)));
            }
          }
        }
        Ltarget[jks] += L;
      }
    }
  }

  /** Kernel L2P operation
   * r += Op(L, t) where L is the local expansion and r is the result
   *
   * @param[in] L The local expansion
   * @param[in] center The center of the box with the local expansion
   * @param[in] target The target of this L2P operation
   * @param[in] result The result to accumulate into
   * @pre L includes the influence of all sources outside its box
   */
  void L2P(const local_type& L, const point_type& center,
           const target_type& target, result_type& result) const {
    complex Ynm[P*P], YnmTheta[P*P];
    point_type r = target - center;
    point_type spherical = point_type();
    point_type cartesian = point_type();
    real rho, theta, phi;
    cart2sph(rho, theta, phi, r);
    evalMultipole(rho, theta, phi, Ynm, YnmTheta);
    for (int n = 0; n < P; ++n) {
      int nm  = n*(n+1);
      int nms = n*(n+1)/2;
      result[0] += std::real(L[nms] * Ynm[nm]);
      spherical[0] += std::real(L[nms] * Ynm[nm]) / rho * n;
      spherical[1] += std::real(L[nms] * YnmTheta[nm]);
      for (int m = 1; m <= n; ++m) {
        nm  = n*(n+1)   + m;
        nms = n*(n+1)/2 + m;
        result[0]    += 2 * std::real(L[nms] * Ynm[nm]);
        spherical[0] += 2 * std::real(L[nms] * Ynm[nm]) / rho * n;
        spherical[1] += 2 * std::real(L[nms] * YnmTheta[nm]);
        spherical[2] += 2 * std::real(L[nms] * Ynm[nm] * complex(0,1)) * m;
      }
    }
    sph2cart(rho, theta, phi, spherical, cartesian);
    result[1] += cartesian[0];
    result[2] += cartesian[1];
    result[3] += cartesian[2];
  }

 protected:

  //! Evaluate solid harmonics \f$ r^n Y_{n}^{m} \f$
  void evalMultipole(real rho, real alpha, real beta,
                     complex *Ynm, complex *YnmTheta) const {
    real x = std::cos(alpha);                                   // x = cos(alpha)
    real y = std::sin(alpha);                                   // y = sin(alpha)
    real fact = 1;                                              // Initialize 2 * m + 1
    real pn = 1;                                                // Initialize Legendre polynomial Pn
    real rhom = 1;                                              // Initialize rho^m
    complex ei = std::exp(complex(0,beta));                     // exp(i * beta)
    complex eim = 1.0;                                          // Initialize exp(i * m * beta)
    for (int m=0; m<P; m++) {                                     // Loop over m in Ynm
      real p = pn;                                              //  Associated Legendre polynomial Pnm
      int npn = m * m + 2 * m;                                    //  Index of Ynm for m > 0
      int nmn = m * m;                                            //  Index of Ynm for m < 0
      Ynm[npn] = rhom * p * eim;                                  //  rho^m * Ynm for m > 0
      Ynm[nmn] = std::conj(Ynm[npn]);                             //  Use conjugate relation for m < 0
      real p1 = p;                                              //  Pnm-1
      p = x * (2 * m + 1) * p1;                                   //  Pnm using recurrence relation
      YnmTheta[npn] = rhom * (p - (m + 1) * x * p1) / y * eim;    // theta derivative of r^n * Ynm
      rhom *= rho;                                                //  rho^m
      real rhon = rhom;                                         //  rho^n
      for (int n=m+1; n<P; n++) {                                 //  Loop over n in Ynm
        int npm = n * n + n + m;                                  //   Index of Ynm for m > 0
        int nmm = n * n + n - m;                                  //   Index of Ynm for m < 0
        rhon /= -(n + m);                                         //   Update factorial
        Ynm[npm] = rhon * p * eim;                                //   rho^n * Ynm
        Ynm[nmm] = std::conj(Ynm[npm]);                           //   Use conjugate relation for m < 0
        real p2 = p1;                                           //   Pnm-2
        p1 = p;                                                   //   Pnm-1
        p = (x * (2 * n + 1) * p1 - (n + m) * p2) / (n - m + 1);  //   Pnm using recurrence relation
        YnmTheta[npm] = rhon * ((n - m + 1) * p - (n + 1) * x * p1) / y * eim;// theta derivative
        rhon *= rho;                                              //   Update rho^n
      }                                                           //  End loop over n in Ynm
      rhom /= -(2 * m + 2) * (2 * m + 1);                         //  Update factorial
      pn = -pn * fact * y;                                        //  Pn
      fact += 2;                                                  //  2 * m + 1
      eim *= ei;                                                  //  Update exp(i * m * beta)
    }                                                             // End loop over m in Ynm
  }

  //! Evaluate singular harmonics \f$ r^{-n-1} Y_n^m \f$
  void evalLocal(real rho, real alpha, real beta,
                 complex *Ynm) const {
    real x = std::cos(alpha);                                   // x = cos(alpha)
    real y = std::sin(alpha);                                   // y = sin(alpha)
    real fact = 1;                                              // Initialize 2 * m + 1
    real pn = 1;                                                // Initialize Legendre polynomial Pn
    real invR = -1.0 / rho;                                     // - 1 / rho
    real rhom = -invR;                                          // Initialize rho^(-m-1)
    complex ei = std::exp(complex(0,beta));                            // exp(i * beta)
    complex eim = 1.0;                                          // Initialize exp(i * m * beta)
    for (int m=0; m<P; m++) {                                     // Loop over m in Ynm
      real p = pn;                                              //  Associated Legendre polynomial Pnm
      int npn = m * m + 2 * m;                                    //  Index of Ynm for m > 0
      int nmn = m * m;                                            //  Index of Ynm for m < 0
      Ynm[npn] = rhom * p * eim;                                  //  rho^(-m-1) * Ynm for m > 0
      Ynm[nmn] = std::conj(Ynm[npn]);                             //  Use conjugate relation for m < 0
      real p1 = p;                                              //  Pnm-1
      p = x * (2 * m + 1) * p1;                                   //  Pnm using recurrence relation
      rhom *= invR;                                               //  rho^(-m-1)
      real rhon = rhom;                                         //  rho^(-n-1)
      for (int n=m+1; n<P; n++) {                                 //  Loop over n in Ynm
        int npm = n * n + n + m;                                  //   Index of Ynm for m > 0
        int nmm = n * n + n - m;                                  //   Index of Ynm for m < 0
        Ynm[npm] = rhon * p * eim;                                //   rho^n * Ynm for m > 0
        Ynm[nmm] = std::conj(Ynm[npm]);                           //   Use conjugate relation for m < 0
        real p2 = p1;                                           //   Pnm-2
        p1 = p;                                                   //   Pnm-1
        p = (x * (2 * n + 1) * p1 - (n + m) * p2) / (n - m + 1);  //   Pnm using recurrence relation
        rhon *= invR * (n - m + 1);                               //   rho^(-n-1)
      }                                                           //  End loop over n in Ynm
      pn = -pn * fact * y;                                        //  Pn
      fact += 2;                                                  //  2 * m + 1
      eim *= ei;                                                  //  Update exp(i * m * beta)
    }                                                             // End loop over m in Ynm
  }

 private:

  /** Spherical to cartesian coordinates */
  inline void sph2cart(real r, real theta, real phi,
                       const point_type& spherical,
                       point_type& cartesian) const {
    real st = std::sin(theta);
    real ct = std::cos(theta);
    real sp = std::sin(phi);
    real cp = std::cos(phi);
    // x component (not x itself)
    cartesian[0] = st * cp * spherical[0]
        + ct * cp / r * spherical[1]
        - sp / r / st * spherical[2];
    // y component (not y itself)
    cartesian[1] = st * sp * spherical[0]
        + ct * sp / r * spherical[1]
        + cp / r / st * spherical[2];
    // z component (not z itself)
    cartesian[2] = ct * spherical[0]
        - st / r * spherical[1];
  }

  /** Cartesian to spherical coordinates */
  inline void cart2sph(real& r, real& theta, real& phi,
                       const point_type& x) const {
    r = norm(x);
    theta = std::acos(x[2] / (r + 1e-100));
    phi = std::atan2(x[1], x[0]);
  }
};
