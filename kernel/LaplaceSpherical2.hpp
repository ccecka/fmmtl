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

#include "kernel/Util/spherical_harmonics.hpp"

class LaplaceSpherical2
    : public fmmtl::Expansion<LaplaceKernel, LaplaceSpherical2> {
 protected:
  typedef double real;
  typedef std::complex<real> complex;

  //! Expansion order
  int P;

  //! (-1)^n
  inline int neg1pow(int n) const {
    return ((n & 1) ? -1 : 1);
  }

 public:
  //! Point type
  typedef Vec<3,real> point_type;

  //! Multipole expansion type
  typedef std::vector<complex> multipole_type;
  //! Local expansion type
  typedef std::vector<complex> local_type;

  //! Constructor
  LaplaceSpherical2(int _P = 5)
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
    real rho, theta, phi;
    cart2sph(rho, theta, phi, r);
    evalMultipole(rho, theta, phi, P, Ynm, YnmTheta);
    for (int n = 0; n < P; ++n) {
      for (int m = 0; m <= n; ++m) {
        int nm  = n*(n+1)   - m;
        int nms = n*(n+1)/2 + m;
        M[nms] += Ynm[nm] * charge;
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
    real rho, theta, phi;
    cart2sph(rho, theta, phi, translation);
    evalMultipole(rho, theta, phi, P, Ynm, YnmTheta);
    for (int j = 0; j < P; ++j) {
      for (int k = 0; k <= j; ++k) {
        int jks = j*(j+1)/2 + k;
        complex M = 0;
        for (int n = 0; n <= j; ++n) {
          for (int m = std::max(-n,-j+k+n); m <= std::min(k-1,n); ++m) {
            int jnkms = (j-n)*(j-n+1)/2 + k - m;
            int nm    = n*(n+1) - m;
            M += Msource[jnkms] * Ynm[nm] * real(neg1pow(m*(m<0)) * neg1pow(n));
          }
          for (int m=k; m<=std::min(n,j+k-n); m++) {
            int jnkms = (j - n) * (j - n + 1) / 2 - k + m;
            int nm    = n * n + n - m;
            M += std::conj(Msource[jnkms]) * Ynm[nm] * real(neg1pow(k+n+m));
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
    complex Ynm[P*P];
    real rho, theta, phi;
    cart2sph(rho, theta, phi, translation);
    evalLocal(rho, theta, phi, P, Ynm);
    for (int j = 0; j < P; ++j) {
      real Cnm = neg1pow(j);
      for (int k = 0; k <= j; k++) {
        int jks = j*(j+1)/2 + k;
        complex L = 0;
        for (int n=0; n<P-j; ++n) {
          for (int m=-n; m<0; ++m) {
            int nms  = n*(n+1)/2 - m;
            int jnkm = (j+n)*(j+n+1) + m - k;
            L += std::conj(Msource[nms]) * Cnm * Ynm[jnkm];
          }
          for (int m = 0; m <= n; ++m) {
            int nms  = n*(n+1)/2 + m;
            int jnkm = (j+n)*(j+n+1) + m - k;
            real Cnm2 = Cnm * neg1pow((k-m)*(k<m)+m);
            L += Msource[nms] * Cnm2 * Ynm[jnkm];
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
    real rho, theta, phi;
    cart2sph(rho, theta, phi, translation);
    evalMultipole(rho, theta, phi, P, Ynm, YnmTheta);
    for (int j = 0; j < P; ++j) {
      for (int k = 0; k <= j; ++k) {
        int jks = j*(j+1)/2 + k;
        complex L = 0;
        for (int n = j; n < P; ++n) {
          for (int m = j+k-n; m < 0; ++m) {
            int jnkm = (n-j)*(n-j+1) + m - k;
            int nms  = n*(n+1)/2 - m;
            L += std::conj(Lsource[nms]) * Ynm[jnkm] * real(neg1pow(k));
          }
          for (int m = 0; m <= n; ++m) {
            if (n-j >= abs(m-k)) {
              int jnkm = (n-j)*(n-j+1) + m - k;
              int nms  = n*(n+1)/2 + m;
              L += Lsource[nms] * Ynm[jnkm] * real(neg1pow((m-k)*(m<k)));
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
    real rho, theta, phi;
    cart2sph(rho, theta, phi, target - center);
    evalMultipole(rho, theta, phi, P, Ynm, YnmTheta);
    point_type spherical = point_type();
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
    point_type cartesian = sph2cart(rho, theta, phi, spherical);
    result[1] += cartesian[0];
    result[2] += cartesian[1];
    result[3] += cartesian[2];
  }
};
