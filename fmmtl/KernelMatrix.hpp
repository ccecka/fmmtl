#pragma once

#include "fmmtl/config.hpp"

//! Global logging
#include "fmmtl/Logger.hpp"
static Logger fmm_global_log;

#include "fmmtl/FMMOptions.hpp"
#include "fmmtl/KernelMatrixPlan.hpp"

#include "fmmtl/meta/kernel_traits.hpp"

#include <vector>

namespace fmmtl {

/** @brief A kernel matrix with kernel of type \c E
 *
 * For a \f$(m \times n)\f$-dimension matrix and
 * \f$ 0 \leq i < m, 0 \leq j < n\f$, every element \f$ m_{i,j} \f$ is generated
 * by the kernel mapping \f$i\f$th target and \f$j\f$th source to a value:
 *
 * @tparam E The expansion type. An expansion is a kernel with additional
 *             operators to accelerate matrix math.
 *
 * TODO: Optimize on aliased source and targets
 */
template <class E>
class kernel_matrix {
  typedef kernel_matrix<E>  this_type;
 public:
  FMMTL_IMPORT_EXPANSION_TRAITS(E);

  typedef std::vector<target_type> target_array;
  typedef std::vector<source_type> source_array;

  typedef std::size_t size_type;

  typedef kernel_value_type  value_type;
  typedef const value_type&  reference;  // ?
  typedef const value_type&  const_reference;

  /** @brief Default constructor. */
  explicit kernel_matrix() : plan(nullptr) {}

  /** @brief         Creates the matrix with the given size
   *
   * @param rows      Number of rows of the matrix
   * @param cols      Number of columns of the matrix
   */
  explicit kernel_matrix(size_type rows, size_type cols)
      : targets_(rows), sources_(cols), plan(nullptr) {}

  template <class TA, class SA>
  explicit kernel_matrix(const expansion_type& e,
                         const TA& targets,
                         const SA& sources)
      : e_(e),
        targets_(targets.begin(), targets.end()),
        sources_(sources.begin(), sources.end()),
        plan(nullptr) {
  }

  template <class SA>
  explicit kernel_matrix(const expansion_type& e,
                         const SA& sources)
      : e_(e),
        targets_(sources.begin(), sources.end()),
        sources_(sources.begin(), sources.end()),
        plan(nullptr) {
  }

  /** @brief Destroy this Kernel matrix */
  ~kernel_matrix() { destroy_plan(); }

  /** @brief Returns the number of rows of the matrix */
  inline size_type rows()  const { return targets_.size(); }
  /** @brief Returns the number of rows of the matrix */
  inline size_type size1() const { return rows(); }

  /** @brief Returns the number of columns of the matrix */
  inline size_type cols()  const { return sources_.size(); }
  /** @brief Returns the number of columns of the matrix */
  inline size_type size2() const { return cols(); }

  /** @brief Returns the kernel generating the elements of the matrix */
  inline       kernel_type& kernel()       { return e_.kernel(); }
  /** @brief Returns the const kernel generating the elements of the matrix */
  inline const kernel_type& kernel() const { return e_.kernel(); }
  /** @brief Returns the expansion operating in the fast MVM */
  inline       expansion_type& expansion()       { return e_.expansion(); }
  /** @brief Returns the const expansion operating in the fast MVM */
  inline const expansion_type& expansion() const { return e_.expansion(); }

  /** @brief Returns a const reference to the source array */
  inline const target_array& targets() const { return targets_; }
  // TODO: Return a proxy?
  /** @brief Returns a const reference to the ith target */
  inline const target_type& target(size_type i) const { return targets_[i]; }
  /** @brief Returns a reference to the ith target */
  inline target_type& target(size_type i) {
    destroy_plan();
    return targets_[i];
  }
  /** @brief Returns a const reference to the permuted target array */
  inline const target_array& permuted_targets() const { return plan->targets(); }

  /** @brief Returns a const reference to the source array */
  inline const source_array& sources() const { return sources_; }
  // TODO: Return a proxy?
  /** @brief Returns a const reference to the jth source */
  inline const source_type& source(size_type j) const { return sources_[j]; }
  /** @brief Returns a reference to the ith target */
  inline source_type& source(size_type j) {
    destroy_plan();
    return sources_[j];
  }
  /** @brief Returns a const reference to the permuted target array */
  inline const target_array& permuted_sources() const { return plan->sources(); }

  /** @brief Returns the matrix element K(i,j) */
  inline value_type operator()(size_type i, size_type j) const {
    FMMTL_ASSERT(i < rows());
    FMMTL_ASSERT(j < cols());
    return expansion()(target(i), source(j));
  }
  /** @brief Syntactic sugar for matrix-vector multiplication
   * TODO: Replace with expression template library */
  inline std::vector<result_type> operator*(const std::vector<charge_type>& c) const {
    std::vector<result_type> result(rows());
    this->prod_impl(c, result);
    return result;
  }

  /** @brief Set options to be passed on to the matrix-vector product plan */
  inline void set_options(const FMMOptions& opts) {
    destroy_plan();
    opts_ = opts;
  }

 private:
  inline void destroy_plan() {
    delete plan;
    plan = nullptr;
  }
  inline void create_plan() const {
    FMMTL_ASSERT(plan == nullptr);
    plan = make_kernel_matrix_plan(*this, opts_);
  }

  template <class CA, class RA>
  inline void prod_impl(const CA& charges, RA& results) const {
    if (plan == nullptr)
      create_plan();
    return plan->execute(charges, results);
  }

  expansion_type e_;
  target_array targets_;
  source_array sources_;

  mutable FMMOptions opts_;
  mutable PlanBase<expansion_type>* plan;
};


} // end namespace fmmtl
