#pragma once

#include "Gaussian.kern"

#include "fmmtl/Expansion.hpp"
#include "fmmtl/numeric/Vec.hpp"
#include "Util/GradedPolynomial.hpp"

// Meh
#include "experimental/util/TensorIndexGridRange.hpp"


/** Gauss expansions in D dimensions using order P polynomials */
template <std::size_t D, std::size_t P>
struct GaussExpansion
    : public fmmtl::Expansion<Gaussian<D>, GaussExpansion<D,P> >
{
  FMMTL_IMPORT_KERNEL_TRAITS(Gaussian<D>);

  typedef Vec<D,double> point_type;

  typedef GradedPolynomial<double,D,P> multipole_type;

  // Not used!
  typedef char local_type;

  // The precomputed Taylor coefficients 2^|alpha| / factorial(alpha)
  std::array<double,multipole_type::size()> coeff;

  // Constructor
  // TODO
  GaussExpansion() {
    // Precompute the Taylor coefficients
    // Dumb way to do this...
    for (auto& i : TensorIndexGridRange<D,P>()) {
      std::size_t idx = poly_index(i);
      if (idx >= coeff.size())
        continue;
      // Compute 2^|alpha| / alpha!
      coeff[idx] = 1;
      for (unsigned k = 0; k < D; ++k)
        coeff[idx] *= std::pow(2.0,i[k]) / factorial(i[k]);
    }
  }

  // XXX: Shouldn't be needed
  void init_multipole(multipole_type& M, const point_type&, unsigned) const {
    M.mono.fill(0);
  }

  // Define the S2M
  // XXX: Should only be leaves
  template <typename SourceIter, typename ChargeIter>
  void S2M(SourceIter s_first, SourceIter s_last, ChargeIter c_first,
           const point_type& center, multipole_type& M) const {
    // For the multi-index alpha,
    // M_alpha = sum_i q_i exp(-R^2/h^2) ((x_i - c) / h)^alpha
    for ( ; s_first != s_last; ++s_first, ++c_first) {
      const source_type& source = *s_first;
      const charge_type& charge = *c_first;
      M += multipole_type(this->inv_h_ * (source - center),
                          charge * this->operator()(source, center));
    }
    // Scale by 2^|alpha| / alpha!
    for (unsigned k = 0; k < M.mono.size(); ++k)
      M.mono[k] *= coeff[k];
  }

  void M2T(const multipole_type& M, const point_type& center,
           const target_type& target, result_type& result) const {
    // q_alpha = exp(-R^2/h^2) ((y_i - c) / h)^alpha
    multipole_type q(this->inv_h_ * (target - center),
                     this->operator()(target, center));
    result += inner_prod(q, M);
  }
};


// TODO: Define custom MAC to select leaf pairs and prune distant pairs
