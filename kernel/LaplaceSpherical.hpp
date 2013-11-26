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
#include "fmmtl/Vec.hpp"

class LaplaceSpherical
    : public fmmtl::Expansion<LaplaceKernel, LaplaceSpherical> {
 protected:
  typedef double real;
  typedef std::complex<real> complex;

  //! Expansion order
  const int P;
  //! \f$ \sqrt{ \frac{(n - |m|)!}{(n + |m|)!} } \f$
  std::vector<real> prefactor;
  //! \f$ (-1)^n / \sqrt{(n + m)! * (n - m)!} \f$
  std::vector<real> Anm;
  //! M2L translation matrix \f$ C_{jn}^{km} \f$
  std::vector<complex> Cnm;

  //! (-1)^n
  inline int neg1pow(int n) const {
    return ((n & 1) ? -1 : 1);
  };
  //! i^n
  inline complex ipow(int n) const {
    switch (n & 3) {
      case 0: return complex(1,0);
      case 1: return complex(0,1);
      case 2: return complex(-1,0);
      case 3: return complex(0,-1);
    }
    assert(false);
    return complex(0,0);
  }
  //! Factorial helper function
  long double factorial(unsigned n) const {
    return std::tgamma(n+1);
  }

  //! Custom multipole type
  struct multipole {
    std::vector<complex> M;
    real RCRIT;
    real RMAX;

    //! Convenience method
    complex& operator[](const int i) {
      return M[i];
    }
    //! Convenience method
    const complex& operator[](const int i) const {
      return M[i];
    }
  };

 public:
  //! The dimension of the spacial interpretation of the source/target_type.
  static const unsigned dimension = 3;
  //! Point type
  typedef Vec<3,real> point_type;

  //! Multipole expansion type
  typedef multipole multipole_type;
  //! Local expansion type
  typedef std::vector<complex> local_type;

  //! Default constructor -- use delegating constructor
  LaplaceSpherical() : LaplaceSpherical(5) {}

  //! Constructor
  LaplaceSpherical(int _P)
      : P(_P), prefactor(4*P*P), Anm(4*P*P), Cnm(P*P*P*P) {
    for (int n = 0; n != 2*P; ++n) {           // Loop over n in Anm
      for (int m = -n; m <= n; ++m) {          //  Loop over m in Anm
        int nm = n*(n+1) + m;                  //   Index of Anm
        // sqrt( (n - |m|)! / (n + |m|)! )
        prefactor[nm] = std::sqrt(factorial(n-abs(m)) / factorial(n+abs(m)));
        // (-1)^n / sqrt( (n + m)! * (n - m)! )
        Anm[nm] = neg1pow(n) / std::sqrt(factorial(n+m) * factorial(n-m));
      }
    }

    for (int j = 0, jk = 0, jknm = 0; j != P; ++j) {
      for (int k = -j; k <= j; ++k, ++jk) {              // jk = j*j + j + m
        for (int n = 0, nm = 0; n != P; ++n) {
          for (int m = -n; m <= n; ++m, ++nm, ++jknm) {  // nm = n*n + n + m?
            int jnkm = (j+n)*(j+n+1) + (m-k);
            Cnm[jknm] = ipow(abs(k-m)-abs(k)-abs(m))
                * (neg1pow(j) * Anm[nm] * Anm[jk] / Anm[jnkm]);
          }
        }
      }
    }
  }

  /** Initialize a multipole expansion with the size of a box at this level */
  void init_multipole(multipole_type& M,
                      const point_type& extents, unsigned) const {
    M.M = std::vector<complex>(P*(P+1)/2, 0);
    M.RMAX = 0;
    M.RCRIT = extents[0] / 2;
  }
  /** Initialize a local expansion with the size of a box at this level */
  void init_local(local_type& L,
                  const point_type&, unsigned) const {
    L = std::vector<complex>(P*(P+1)/2, 0);
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
    complex Ynm[4*P*P], YnmTheta[4*P*P];
    point_type dist = source - center;
    real rho, alpha, beta;
    cart2sph(rho,alpha,beta,dist);
    evalMultipole(rho,alpha,-beta,Ynm,YnmTheta);
    for (int n = 0; n != P; ++n) {
      for (int m = 0; m <= n; ++m) {
        const int nm  = n*(n+1)   + m;
        const int nms = n*(n+1)/2 + m;
        M[nms] += charge * Ynm[nm];
      }
    }
    M.RMAX = std::max(M.RMAX, norm(dist));
    M.RCRIT = std::min(M.RCRIT, M.RMAX);
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
    complex Ynm[4*P*P], YnmTheta[4*P*P];
    real Rmax = Mtarget.RMAX;
    real R = norm(translation) + Msource.RCRIT;
    if (R > Rmax) Rmax = R;
    real rho, alpha, beta;
    cart2sph(rho,alpha,beta,translation);
    evalMultipole(rho,alpha,-beta,Ynm,YnmTheta);
    for( int j=0; j!=P; ++j ) {
      for( int k=0; k<=j; ++k ) {
        const int jk = j*(j+1)    + k;
        const int jks = j*(j+1)/2 + k;
        complex M = 0;
        for( int n=0; n<=j; ++n ) {
          for( int m=-n; m<=std::min(k-1,n); ++m ) {
            if( j-n >= k-m ) {
              const int jnkm  = (j-n)*(j-n+1)   + (k-m);
              const int jnkms = (j-n)*(j-n+1)/2 + (k-m);
              const int nm    = n*(n+1) + m;
              M += Msource[jnkms] * ipow(m-abs(m)) * Ynm[nm]
              * real(neg1pow(n) * Anm[nm] * Anm[jnkm] / Anm[jk]);
            }
          }
          for( int m=k; m<=n; ++m ) {
            if( j-n >= m-k ) {
              const int jnkm  = (j-n)*(j-n+1)   + (k-m);
              const int jnkms = (j-n)*(j-n+1)/2 - (k-m);
              const int nm    = n*(n+1) + m;
              M += std::conj(Msource[jnkms]) * Ynm[nm]
              * real(neg1pow(k+n+m) * Anm[nm] * Anm[jnkm] / Anm[jk]);
            }
          }
        }
        Mtarget[jks] += M;
      }
    }
    Mtarget.RMAX = Rmax;
    Mtarget.RCRIT = std::min(Mtarget.RCRIT, Mtarget.RMAX);
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
    complex Ynm[4*P*P], YnmTheta[4*P*P];

    point_type dist = translation;
    real rho, alpha, beta;
    cart2sph(rho,alpha,beta,dist);
    evalLocal(rho,alpha,beta,Ynm,YnmTheta);
    for( int j=0; j!=P; ++j ) {
      for( int k=0; k<=j; ++k ) {
        const int jk = j * j + j + k;
        const int jks = j * (j + 1) / 2 + k;
        complex L = 0;
        for( int n=0; n!=P; ++n ) {
          for( int m=-n; m<0; ++m ) {
            const int nm   = n*(n+1)   + m;
            const int nms  = n*(n+1)/2 - m;
            const int jknm = jk*P*P + nm;
            const int jnkm = (j+n)*(j+n+1) + (m-k);
            L += std::conj(Msource[nms]) * Cnm[jknm] * Ynm[jnkm];
          }
          for( int m=0; m<=n; ++m ) {
            const int nm   = n*(n+1)   + m;
            const int nms  = n*(n+1)/2 + m;
            const int jknm = jk*P*P + nm;
            const int jnkm = (j+n)*(j+n+1) + (m-k);
            L += Msource[nms] * Cnm[jknm] * Ynm[jnkm];
          }
        }
        Ltarget[jks] += L;
      }
    }
  }

  /** Kernel M2P operation
   * r += Op(M, t) where M is the multipole and r is the result
   *
   * @param[in] M The multpole expansion
   * @param[in] center The center of the box with the multipole expansion
   * @param[in] target The target to evaluate the multipole at
   * @param[in,out] result The target's corresponding result to accumulate into
   * @pre M includes the influence of all sources within its box
   */
  void M2P(const multipole_type& M, const point_type& center,
           const target_type& target, result_type& result) const {
    complex Ynm[4*P*P], YnmTheta[4*P*P];
    point_type dist = target - center;
    point_type spherical = point_type();
    point_type cartesian = point_type();
    real r, theta, phi;
    cart2sph(r,theta,phi,dist);
    evalLocal(r,theta,phi,Ynm,YnmTheta);
    for( int n=0; n!=P; ++n ) {
      int nm  = n*(n+1);
      int nms = n*(n+1)/2;
      result[0] += std::real(M[nms] * Ynm[nm]);
      spherical[0] -= std::real(M[nms] * Ynm[nm]) / r * (n+1);
      spherical[1] += std::real(M[nms] * YnmTheta[nm]);
      for( int m=1; m<=n; ++m ) {
        nm  = n*(n+1)   + m;
        nms = n*(n+1)/2 + m;
        result[0] += 2 * std::real(M[nms] * Ynm[nm]);
        spherical[0] -= 2 * std::real(M[nms] *Ynm[nm]) / r * (n+1);
        spherical[1] += 2 * std::real(M[nms] *YnmTheta[nm]);
        spherical[2] += 2 * std::real(M[nms] *Ynm[nm] * complex(0,1)) * m;
      }
    }
    sph2cart(r,theta,phi,spherical,cartesian);
    result[1] += cartesian[0];
    result[2] += cartesian[1];
    result[3] += cartesian[2];
  }

  /** Kernel L2L operator
   * L_t += Op(L_s) where L_t is the target and L_s is the source
   *
   * @param[in] source The local source at the parent level
   * @param[in,out] target The local target to accumulate into
   * @param[in] translation The vector from source to target
   * @pre Lsource includes the influence of all points outside its box
   */
  void L2L(const local_type& source,
           local_type& target,
           const point_type& translation) const {
    complex Ynm[4*P*P], YnmTheta[4*P*P];
    real rho, alpha, beta;
    cart2sph(rho,alpha,beta,translation);
    evalMultipole(rho,alpha,beta,Ynm,YnmTheta);
    for( int j=0; j!=P; ++j ) {
      for( int k=0; k<=j; ++k ) {
        const int jk  = j*(j+1)   + k;
        const int jks = j*(j+1)/2 + k;
        complex L = 0;
        for( int n=j; n!=P; ++n ) {
          for( int m=j+k-n; m<0; ++m ) {
            const int jnkm = (n-j)*(n-j+1) + (m-k);
            const int nm   = n*(n+1)   - m;
            const int nms  = n*(n+1)/2 - m;
            L += std::conj(source[nms]) * Ynm[jnkm]
                * real(neg1pow(k) * Anm[jnkm] * Anm[jk] / Anm[nm]);
          }
          for( int m=0; m<=n; ++m ) {
            if( n-j >= abs(m-k) ) {
              const int jnkm = (n-j)*(n-j+1) + (m-k);
              const int nm   = n*(n+1)   + m;
              const int nms  = n*(n+1)/2 + m;
              L += source[nms] * ipow((m-k)-abs(m-k))
                  * Ynm[jnkm] * Anm[jnkm] * Anm[jk] / Anm[nm];
            }
          }
        }
        target[jks] += L;
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
    complex Ynm[4*P*P], YnmTheta[4*P*P];
    point_type dist = target - center;
    point_type spherical = point_type();
    point_type cartesian = point_type();
    real r, theta, phi;
    cart2sph(r,theta,phi,dist);
    evalMultipole(r,theta,phi,Ynm,YnmTheta);
    for( int n=0; n!=P; ++n ) {
      int nm  = n*(n+1);
      int nms = n*(n+1)/2;
      result[0] += std::real(L[nms] * Ynm[nm]);
      spherical[0] += std::real(L[nms] * Ynm[nm]) / r * n;
      spherical[1] += std::real(L[nms] * YnmTheta[nm]);
      for( int m=1; m<=n; ++m ) {
        nm  = n*(n+1)   + m;
        nms = n*(n+1)/2 + m;
        result[0] += 2 * std::real(L[nms] * Ynm[nm]);
        spherical[0] += 2 * std::real(L[nms] * Ynm[nm]) / r * n;
        spherical[1] += 2 * std::real(L[nms] * YnmTheta[nm]);
        spherical[2] += 2 * std::real(L[nms] * Ynm[nm] * complex(0,1)) * m;
      }
    }
    sph2cart(r,theta,phi,spherical,cartesian);
    result[1] += cartesian[0];
    result[2] += cartesian[1];
    result[3] += cartesian[2];
  }

 protected:

  //! Evaluate solid harmonics \f$ r^n Y_{n}^{m} \f$
  void evalMultipole(real rho, real alpha, real beta,
                     complex* Ynm, complex* YnmTheta) const {
    return generate_Y(rho, alpha, beta,
                      Ynm, YnmTheta,
                      1, P);
  }

  //! Generate singular harmonics \f$ r^{-n-1} Y_n^m \f$
  void evalLocal(real rho, real alpha, real beta,
                 complex* Ynm, complex* YnmTheta) const {
    return generate_Y(1.0/rho, alpha, beta,
                      Ynm, YnmTheta,
                      1.0/rho, 2*P);
  }

 private:

  //! Generate singular harmonics \f$ A * r^{n} Y_n^m \f$ and theta derivative
  inline void generate_Y(real rho, real alpha, real beta,
                         complex* Ynm, complex* YnmTheta,
                         real A, int K) const {
    real x = std::cos(alpha);                        // x = cos(alpha)
    real y = std::sin(alpha);                        // y = sin(alpha)
    real pn = 1;                                     // Initialize Legendre Pn
    for (int m = 0; m != K; ++m) {                   // Loop over m in Ynm
      complex eimb = std::exp(complex(0,m*beta));    //  exp(i * m * beta)
      int npn = m*(m+2);                             //  Index of Ynm for m > 0
      int nmn = m*m;                                 //  Index of Ynm for m < 0
      Ynm[npn] = A * pn * prefactor[npn] * eimb;     //  rho^m * Ynm for m > 0
      Ynm[nmn] = std::conj(Ynm[npn]);                //  Use conj for m < 0
      YnmTheta[npn] = m*x/y * Ynm[npn];              // theta derivative

      real p1 = pn;                                  //  Pnm-1
      real p = (2*m+1) * x * pn;                     //  Pnm using recurrence
      real rhon = A;                                 //  rho^n
      for (int n = m+1; n != K; ++n) {               //  Loop over n in Ynm
        rhon *= rho;                                 //   Update rho^n
        int npm = n*(n+1) + m;                       //   Index of Ynm for m > 0
        int nmm = n*(n+1) - m;                       //   Index of Ynm for m < 0
        Ynm[npm] = rhon * p * prefactor[npm] * eimb; //   rho^n * Ynm
        Ynm[nmm] = std::conj(Ynm[npm]);              //   Use conj for m < 0
        YnmTheta[npm] = (n*x - (n+m)*p1/p) / y * Ynm[npm];  // theta derivative
        real p2 = p1;                                //   Pnm-2
        p1 = p;                                      //   Pnm-1
        p = ((2*n+1) * x * p1 - (n+m) * p2) / (n-m+1); //   Pnm using recurrence
      }                                              //  End loop over n in Ynm
      A *= rho;                                      //  rho^m
      pn = -pn * (2*m+1) * y;                        //  Pn
    }                                                // End loop over m in Ynm
  }

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
    r = std::sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]);
    theta = std::acos(x[2] / (r + 1e-100));
    phi = std::atan2(x[1], x[0]);
  }
};
