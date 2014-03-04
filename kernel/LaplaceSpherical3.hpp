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

class LaplaceSpherical3
    : public fmmtl::Expansion<LaplaceKernel, LaplaceSpherical3> {
 protected:
  typedef double real;
  typedef std::complex<real> complex;

  //! Expansion order
  int P;

  //! (-1)^n
  inline static constexpr real neg1pow(int n) {
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
  LaplaceSpherical3(int _P = 5)
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
    complex Ynm[P*P];
    real rho, theta, phi;
    cart2sph(rho, theta, phi, center - source);
    evalZ(rho, theta, phi, P, Ynm);
    for (int n = 0; n < P; ++n) {
      for (int m = 0; m <= n; ++m) {
        int nm  = n*(n+1)/2 + m;
        M[nm] += neg1pow(m) * std::conj(Ynm[nm]) * charge;
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
    complex Y[P*P];
    real rho, theta, phi;
    cart2sph(rho, theta, phi, translation);
    evalZ(rho, theta, phi, P, Y);
    for (int n = 0; n != P; ++n) {
      for (int m = 0; m <= n; ++m) {
        complex M = 0;
        for (int j = 0; j <= n; ++j) {
          // All k with -j <= k <= 0 and 0 <= m-k <= n-j
          // Thus, k >= -j and k >= -n+j+m
          int k = std::max(-j, -n+j+m);
          // Thus, k <= 0 and k <= m
          for ( ; k <= 0; ++k) {
            // k is negative and m-k is positive
            int midx = (n-j)*(n-j+1)/2 + (m-k);
            int jk   = j*(j+1)/2 - k;
            M += Y[jk] * Msource[midx];
          }

          // All k with 0 < k <= j and 0 <= m-k <= n-j
          // Thus, k <= j and k <= m
          int end = std::min(j, m);
          for ( ; k <= end; ++k) {
            // k is positive and m-k is positive
            int midx = (n-j)*(n-j+1)/2 + (m-k);
            int jk   = j*(j+1)/2 + k;
            M += neg1pow(k) * std::conj(Y[jk]) * Msource[midx];
          }

          // All k with 0 <= k < j and -(n-j) <= m-k <= 0
          // Thus, k <= j and k <= n-j+m
          end = std::min(j, n-j+m);
          for ( ; k <= end; ++k) {
            // k is positive and m-k is negative
            int midx = (n-j)*(n-j+1)/2 - (m-k);
            int jk   = j*(j+1)/2 + k;
            M += std::conj(neg1pow(m) * Y[jk] * Msource[midx]);
          }
        }
        int nm = n*(n+1)/2 + m;
        Mtarget[nm] += M;
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
    complex Y[4*P*P];
    real rho, theta, phi;
    cart2sph(rho, theta, phi, translation);
    evalW(rho, theta, phi, 2*P, Y);
    for (int n = 0; n != P; ++n) {
      for (int m = 0; m <= n; ++m) {
        complex L = 0;
        for (int j = 0; j != P; ++j) {
          // All k with -j <= k <= 0 and -(j+n) <= k-m <= 0
          // Thus, k >= -j and k >= m-n-j
          int k = -j;
          // Thus, k <= 0 and k <= m
          for ( ; k <= 0; ++k) {
            // k is negative and k-m is negative
            int midx = j*(j+1)/2 - k;
            int yidx = (j+n)*(j+n+1)/2 - (k-m);
            L += std::conj(neg1pow(m) * Msource[midx] * Y[yidx]);
          }

          // All k with 0 <= k <= j and -(j+n) <= k-m <= 0
          // Thus, k <= j and k <= m
          int end = std::min(j, m);
          for ( ; k <= end; ++k) {
            // k is positive and k-m is negative
            int midx = j*(j+1)/2 + k;
            int yidx = (j+n)*(j+n+1)/2 - (k-m);
            L += neg1pow(k-m) * Msource[midx] * std::conj(Y[yidx]);
          }

          // All k with 0 <= k <= j and 0 <= k-m <= j+n
          // Thus, k <= j and k <= m+n+j
          for ( ; k <= j; ++k) {
            // k is positive and k-m is positive
            int midx = j*(j+1)/2 + k;
            int yidx = (j+n)*(j+n+1)/2 + (k-m);
            L += Msource[midx] * Y[yidx];
          }
        }
        int nm = n*(n+1)/2 + m;
        Ltarget[nm] += L;
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
    complex Y[P*P];
    real rho, theta, phi;
    cart2sph(rho, theta, phi, translation);
    evalZ(rho, theta, phi, P, Y);
    for (int n = 0; n != P; ++n) {
      for (int m = 0; m <= n; ++m) {
        complex L = 0;
        for (int j = n; j != P; ++j) {
          // All k with -j <= k <= 0 and -(j-n) <= k-m <= 0
          // Thus, k >= -j and k >= n+m-j
          int k = n+m-j;
          // Thus, k <= 0 and k <= m
          for ( ; k <= 0; ++k) {
            // k is negative and k-m is negative
            int lidx = j*(j+1)/2 - k;
            int yidx = (j-n)*(j-n+1)/2 - (k-m);
            L += std::conj(neg1pow(m) * Lsource[lidx] * Y[yidx]);
          }

          // All k with 0 <= k <= j and -(j-n) <= k-m <= 0
          // Thus, k <= j and k <= m
          int end = std::min(j, m);
          for ( ; k <= end; ++k) {
            // k is positive and k-m is negative
            int lidx = j*(j+1)/2 + k;
            int yidx = (j-n)*(j-n+1)/2 - (k-m);
            L += neg1pow(k-m) * Lsource[lidx] * std::conj(Y[yidx]);
          }

          // All k with 0 <= k <= j and 0 <= k-m <= j-n
          // Thus, k <= j and k <= m-n+j
          end = std::min(j, m-n+j);
          for ( ; k <= end; ++k) {
            // k is positive and k-m is positive
            int lidx = j*(j+1)/2 + k;
            int yidx = (j-n)*(j-n+1)/2 + (k-m);
            L += Lsource[lidx] * Y[yidx];
          }
        }
        int nm = n*(n+1)/2 + m;
        Ltarget[nm] += L;
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
    evalZ(rho, theta, phi, P, Ynm, YnmTheta);
    point_type spherical = point_type();
    for (int n = 0; n != P; ++n) {
      int nm  = n*(n+1)/2;
      result[0]    += std::real(L[nm] * Ynm[nm]);
      spherical[0] += std::real(L[nm] * Ynm[nm]) / rho * n;
      spherical[1] += std::real(L[nm] * YnmTheta[nm]);
      for (int m = 1; m <= n; ++m) {
        nm  = n*(n+1)/2 + m;
        result[0]    += 2 * std::real(L[nm] * Ynm[nm]);
        spherical[0] += 2 * std::real(L[nm] * Ynm[nm]) / rho * n;
        spherical[1] += 2 * std::real(L[nm] * YnmTheta[nm]);
        spherical[2] += 2 * std::real(L[nm] * Ynm[nm] * complex(0,1)) * m;
      }
    }
    point_type cartesian = sph2cart(rho, theta, phi, spherical);
    result[1] += cartesian[0];
    result[2] += cartesian[1];
    result[3] += cartesian[2];
  }
};
