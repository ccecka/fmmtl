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
  typedef double real_type;
  typedef std::complex<real_type> complex_type;

  //! Expansion order
  int P;

  //! (-1)^n
  inline static constexpr real_type neg1pow(int n) {
    return ((n & 1) ? -1 : 1);
  }

 public:
  //! Point type
  typedef Vec<3,real_type> point_type;

  //! Multipole expansion type
  typedef std::vector<complex_type> multipole_type;
  //! Local expansion type
  typedef std::vector<complex_type> local_type;

  //! Constructor
  LaplaceSpherical3(int _P = 5)
      : P(_P) {
  }

  /** Initialize a multipole expansion with the size of a box at this level */
  void init_multipole(multipole_type& M,
                      const point_type&, unsigned) const {
    M = std::vector<complex_type>(P*(P+1)/2);
  }
  /** Initialize a local expansion with the size of a box at this level */
  void init_local(local_type& L,
                  const point_type&, unsigned) const {
    L = std::vector<complex_type>(P*(P+1)/2);
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
    real_type rho, theta, phi;
    cart2sph(rho, theta, phi, center - source);
    complex_type Z[P*(P+1)/2];   // Avoid initialization?
    evalZ(rho, theta, phi, P, Z);
    int nm = 0;   // n*(n+1)/2+m
    for (int n = 0; n < P; ++n) {
      for (int m = 0; m <= n; ++m, ++nm) {
        M[nm] += neg1pow(m) * std::conj(Z[nm]) * charge;
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
    real_type rho, theta, phi;
    cart2sph(rho, theta, phi, translation);
    complex_type Z[P*(P+1)/2];
    evalZ(rho, theta, phi, P, Z);
    int nm = 0;   // n*(n+1)/2+m
    for (int n = 0; n != P; ++n) {
      for (int m = 0; m <= n; ++m, ++nm) {
        complex_type& M = Mtarget[nm];

        for (int j = 0; j <= n; ++j) {
          // Compute the offset for Y_j
          const auto Zj = Z + j*(j+1)/2;
          // Compute the offset for M_{n-j}
          const auto Mnj = Msource.begin() + (n-j)*(n-j+1)/2;

          // All k with -j <= k <= 0 and 0 <= m-k <= n-j
          // Thus, k >= -j and k >= -n+j+m
          int k = std::max(-j, -n+j+m);
          // Thus, k <= 0 and k <= m
          for ( ; k <= 0; ++k) {
            // k is negative and m-k is positive
            M += Zj[-k] * Mnj[m-k];
          }

          // All k with 0 < k <= j and 0 <= m-k <= n-j
          // Thus, k <= j and k <= m
          int end = std::min(j, m);
          for ( ; k <= end; ++k) {
            // k is positive and m-k is positive
            M += neg1pow(k) * std::conj(Zj[k]) * Mnj[m-k];
          }

          // All k with 0 <= k < j and -(n-j) <= m-k <= 0
          // Thus, k <= j and k <= n+m-j
          end = std::min(j, n+m-j);
          for ( ; k <= end; ++k) {
            // k is positive and m-k is negative
            M += std::conj(neg1pow(m) * Zj[k] * Mnj[k-m]);
          }
        }
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
    real_type rho, theta, phi;
    cart2sph(rho, theta, phi, translation);
    complex_type W[P*(2*P+1)];
    evalW(rho, theta, phi, 2*P, W);
    int nm = 0;    // n*(n+1)/2 + m
    for (int n = 0; n != P; ++n) {
      for (int m = 0; m <= n; ++m, ++nm) {
        complex_type& L = Ltarget[nm];

        for (int j = 0; j != P; ++j) {
          // Compute the offset for M_j
          const auto Mj = Msource.begin() + j*(j+1)/2;
          // Compute the offset for W_{j+n}
          const auto Wjn = W + (j+n)*(j+n+1)/2;

          // All k with -j <= k <= 0 and -(j+n) <= k-m <= 0
          // Thus, k >= -j and k >= m-n-j
          int k = -j;
          // Thus, k <= 0 and k <= m
          for ( ; k <= 0; ++k) {
            // k is negative and k-m is negative
            L += std::conj(neg1pow(m) * Mj[-k] * Wjn[m-k]);
          }

          // All k with 0 <= k <= j and -(j+n) <= k-m <= 0
          // Thus, k <= j and k <= m
          int end = std::min(j, m);
          for ( ; k <= end; ++k) {
            // k is positive and k-m is negative
            L += neg1pow(k-m) * Mj[k] * std::conj(Wjn[m-k]);
          }

          // All k with 0 <= k <= j and 0 <= k-m <= j+n
          // Thus, k <= j and k <= m+n+j
          for ( ; k <= j; ++k) {
            // k is positive and k-m is positive
            L += Mj[k] * Wjn[k-m];
          }
        }
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
    real_type rho, theta, phi;
    cart2sph(rho, theta, phi, translation);
    complex_type Z[P*(P+1)/2];
    evalZ(rho, theta, phi, P, Z);
    int nm = 0;    // n*(n+1)/2 + m
    for (int n = 0; n != P; ++n) {
      for (int m = 0; m <= n; ++m, ++nm) {
        complex_type& L = Ltarget[nm];

        for (int j = n; j != P; ++j) {
          // Compute the offset for L_j
          const auto Lj = Lsource.begin() + j*(j+1)/2;
          // Compute the offset for Z_{j-n}
          const auto Zjn = Z + (j-n)*(j-n+1)/2;

          // All k with -j <= k <= 0 and -(j-n) <= k-m <= 0
          // Thus, k >= -j and k >= n+m-j
          int k = n+m-j;
          // Thus, k <= 0 and k <= m
          for ( ; k <= 0; ++k) {
            // k is negative and k-m is negative
            L += std::conj(neg1pow(m) * Lj[-k] * Zjn[m-k]);
          }

          // All k with 0 <= k <= j and -(j-n) <= k-m <= 0
          // Thus, k <= j and k <= m
          int end = std::min(j, m);
          for ( ; k <= end; ++k) {
            // k is positive and k-m is negative
            L += neg1pow(k-m) * Lj[k] * std::conj(Zjn[m-k]);
          }

          // All k with 0 <= k <= j and 0 <= k-m <= j-n
          // Thus, k <= j and k <= m-n+j
          end = std::min(j, m-n+j);
          for ( ; k <= end; ++k) {
            // k is positive and k-m is positive
            L += Lj[k] * Zjn[k-m];
          }
        }
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
    using std::real;
    using std::imag;

    real_type rho, theta, phi;
    cart2sph(rho, theta, phi, target - center);
    complex_type Z[P*(P+1)/2], dZ[P*(P+1)/2];
    evalZ(rho, theta, phi, P, Z, dZ);

    point_type sph = point_type();
    int nm = 0;
    for (int n = 0; n != P; ++n) {
      const real_type LZ = real(L[nm])*real(Z[nm]) - imag(L[nm])*imag(Z[nm]);
      result[0] += LZ;
      sph[0]    += LZ / rho * n;
      sph[1]    += real(L[nm])*real(dZ[nm]) - imag(L[nm])*imag(dZ[nm]);

      ++nm;
      for (int m = 1; m <= n; ++m, ++nm) {
        const complex_type LZ = L[nm] * Z[nm];
        result[0] += 2 * std::real(LZ);
        sph[0]    += 2 * std::real(LZ) / rho * n;
        sph[1]    += 2 * (real(L[nm])*real(dZ[nm])-imag(L[nm])*imag(dZ[nm]));
        sph[2]    += 2 *-std::imag(LZ) * m;
      }
    }

    const point_type cart = sph2cart(rho, theta, phi, sph);
    result[1] += cart[0];
    result[2] += cart[1];
    result[3] += cart[2];
  }
};
