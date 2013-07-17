#pragma once
/** @file Helmholtz
 * @brief An implementation of the Helmholtz kernel for FMM/treecodes
 *
 * This class will be used to evaluate the matrix-vector product
 * r_i = sum_j K(t_i, s_j) c_j
 *
 * where K is the Helmholtz Green's function:
 * K(t,s) = exp(i kappa |s-t|) / |s-t|
 *
 * and the source_type and target_type are 3D points,
 * and the charge_type, result_type, and kernel_value_type are complex.
 */

#include "fmmtl/Kernel.hpp"

#include "Helmholtz.kern"

#include "Util/HelmholtzQuad.hpp"
#include "Util/HelmholtzReterp.hpp"

#include <complex>
#include <memory>  // using std::shared_ptr

class HelmholtzFourier
    : public HelmholtzPotential,
      public FMM_Expansion<HelmholtzFourier>
{
  typedef double real;
  typedef std::complex<real> complex;

  template <typename T>
  inline complex unit_polar(T phi) const {
    return complex(std::cos(phi), std::sin(phi));
  }

 public:
  //! Point type -- 3D point
  typedef Vec<3,double> point_type;

  //! Multipole expansion type
  typedef S2Function multipole_type;
  //! Local expansion type
  typedef S2Function local_type;

  //! Error bound for the Helmholtz kernel
  double eps;
  //! Algorithmic parameter for accuracy looseness
  double alpha;

  //! Quadrature for each level
  mutable std::vector<std::shared_ptr<S2Quadrature>> level2quad;
  //! Array of reterpolators from level to level-1
  mutable std::vector<std::shared_ptr<Reterpolator>> upReterp;
  //! Array of reterpolators from level to level+1
  mutable std::vector<std::shared_ptr<Reterpolator>> downReterp;

  /** Constructor */
  Helmholtz(double _kappa, double _eps = 1e-4, double _alpha = 0.9)
      : kappa(_kappa), eps(_eps), alpha(_alpha),
        level2quad(20,nullptr), upReterp(20,nullptr), downReterp(20,nullptr) {
  }
  /** Initialize a multipole expansion with the size of a box and level number
   *
   * @param[in] M The multipole to be initialized
   * @param[in] extents The dimensions of the box containing the multipole
   * @param[in] level The level number of the box. 0: Root box
   */
  void init_multipole(multipole_type& M,
                      const point_type& extents, unsigned level) const {
    if (M.quad == nullptr) {
      if (level2quad[level] == nullptr) {
        // Construct the quadrature for this level
        level2quad[level]
            = std::make_shared<S2Quadrature>(kappa, level, extents[0],
                                             alpha, eps);
        std::cerr << *level2quad[level] << std::endl;
        // Check the level above
        if (upReterp[level] == nullptr && level2quad[level-1] != nullptr) {
          // Construct the up interpolator
          upReterp[level]
              = std::make_shared<Reterpolator>(level2quad[level].get(),
                                               level2quad[level-1].get());
          // Construct the down interpolator
          downReterp[level-1]
              = std::make_shared<Reterpolator>(level2quad[level-1].get(),
                                               level2quad[level].get());
        }
        // Check the level below
        if (downReterp[level] == nullptr && level2quad[level+1] != nullptr) {
          // Construct the down reterpolator
          downReterp[level]
              = std::make_shared<Reterpolator>(level2quad[level].get(),
                                               level2quad[level+1].get());
          upReterp[level+1]
              = std::make_shared<Reterpolator>(level2quad[level+1].get(),
                                               level2quad[level].get());
        }
      }
      M = multipole_type(level2quad[level].get());
    }
  }
  /** Initialize a local expansion with the size of a box at this level
   *
   * @param[in] M The multipole to be initialized
   * @param[in] extents The dimensions of the box containing the multipole
   * @param[in] level The level number of the box. 0: Root box
   */
  void init_local(local_type& L,
                  const point_type& extents, unsigned level) const {
    init_multipole(L, extents, level);
  }

  /** Kernel P2M operation
   * M += Op(s) * c where M is the multipole and s is the source
   *
   * @param[in] source The source to accumulate into the multipole
   * @param[in] charge The source's corresponding charge
   * @param[in] center The center of the box containing the multipole expansion
   * @param[in,out] M The multipole expansion to accumulate into
   */
  void P2M(const source_type& source, const charge_type& charge,
           const point_type& center, multipole_type& M) const {
    point_type r = center - source;

    const S2Quadrature& quad = *(M.quad);
    unsigned K = M.size();
    for (unsigned k = 0; k < K; ++k)
      M[k] += charge * unit_polar(kappa * r.dot(quad[k]));
  }

  /** Kernel M2M operator
   * M_t += Op(M_s) where M_t is the target and M_s is the source
   *
   * @param[in] source The multipole source at the child level
   * @param[in,out] target The multipole target to accumulate into
   * @param[in] translation The vector from source to target
   * @pre Msource includes the influence of all sources within its box
   */
  void M2M(const multipole_type& source,
           multipole_type& target,
           const point_type& translation) const {
    assert(source.quad != nullptr);
    assert(target.quad != nullptr);
    int level = source.quad->level;
    assert(target.quad->level == level-1);
    assert(upReterp[level] != nullptr);

    // TODO: optimize
    multipole_type temp = target;
    upReterp[level]->apply(source, temp);

    const S2Quadrature& quad = *(target.quad);
    unsigned K = target.size();
    for (unsigned k = 0; k < K; ++k)
      target[k] += temp[k] * unit_polar(kappa * translation.dot(quad[k]));
  }

  /** Kernel M2L operation
   * L += Op(M) where L is the local expansion and M is the multipole
   *
   * @param[in] source The multpole expansion source
   * @param[in,out] target The local expansion target
   * @param[in] translation The vector from source to target
   * @pre translation obeys the multipole-acceptance criteria
   * @pre source includes the influence of all sources within its box
   */
  void M2L(const multipole_type& source,
           local_type& target,
           const point_type& translation) const {
    assert(source.quad == target.quad);

    // TODO: Precompute by translation vector and Optimize
    Transfer_Function T(source.quad, translation);

    unsigned K = target.size();
    for (unsigned k = 0; k < K; ++k)
      target[k] += T[k] * source[k];
  }

  // XXX! TESTING ONLY
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
    assert(M.quad != nullptr);
    point_type r = target - center;
    Transfer_Function T(M.quad, r);

    result_type I = complex(0,0);
    const S2Quadrature& quad = *(M.quad);
    unsigned K = M.size();
    for (unsigned k = 0; k < K; ++k)
      I += quad.weight(k) * M[k] * T[k];

    result += I;
  }

  /** Kernel L2L operator
   * L_t += Op(L_s) where L_t is the target and L_s is the source
   *
   * @param[in] source The local source at the parent level
   * @param[in,out] target The local target to accumulate into
   * @param[in] translation The vector from source to target
   * @pre source includes the influence of all sources outside its box
   */
  void L2L(const local_type& source,
           local_type& target,
           const point_type& translation) const {
    assert(source.quad != nullptr);
    assert(target.quad != nullptr);
    int level = source.quad->level;
    assert(target.quad->level == level+1);
    assert(downReterp[level] != nullptr);

    // TODO: optimize
    local_type temp = target;
    downReterp[level]->apply(source, temp);

    const S2Quadrature& quad = *(target.quad);
    unsigned K = target.size();
    for (unsigned k = 0; k < K; ++k)
      target[k] += temp[k] * unit_polar(kappa * translation.dot(quad[k]));
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
    point_type r = target - center;
    result_type I = complex(0,0);

    const S2Quadrature& quad = *(L.quad);
    unsigned K = L.size();
    for (unsigned k = 0; k < K; ++k)
      I += quad.weight(k) * L[k] * unit_polar(kappa * r.dot(quad[k]));

    result += I;
  }
};
