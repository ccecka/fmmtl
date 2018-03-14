#pragma once

#include <type_traits>
#include <iterator>

#include <thrust/iterator/transform_iterator.h>
#include <thrust/detail/raw_reference_cast.h>

#include "fmmtl/meta/func_traits.hpp"


// The thrust::transform_iterator does not provide a mutable iterator
// This prevents things like a transform_iterator<identity, > from acting as
// expected and potentially creating uncessary copies.

// Define a thrust::v2::transform_iterator that avoids these copies.
// The caveat is that this iterator may only be dereferenced in the space the
// underlying iterator traverses.

namespace thrust {

namespace v2 {

// See thrust::transform_iterator
template <class AdaptableUnaryFunction,
          class Iterator,
          class Reference = decltype(std::declval<AdaptableUnaryFunction>()(*std::declval<Iterator>())),
          class Value = typename std::remove_cv<Reference>::type>
class transform_iterator
    : public thrust::iterator_adaptor<transform_iterator<AdaptableUnaryFunction,
                                                         Iterator,
                                                         Reference,
                                                         Value>,
                                      Iterator,
                                      Value,
                                      thrust::use_default,   // Leave the system alone
                                      typename thrust::iterator_traits<Iterator>::iterator_category,
                                      Reference>
{
  /*! \cond
   */
 public:
  typedef typename thrust::iterator_adaptor<transform_iterator,
                                            Iterator,
                                            Value,
                                            thrust::use_default,   // Leave the system alone
                                            typename thrust::iterator_traits<Iterator>::iterator_category,
                                            Reference>
  super_t;

  friend class thrust::iterator_core_access;
  /*! \endcond
   */

 public:
  /*! Null constructor does nothing.
   */
  __host__ __device__
  transform_iterator() {}

  /*! This constructor takes as arguments an \c Iterator and an \c AdaptableUnaryFunction
   *  and copies them to a new \p transform_iterator.
   *
   *  \param x An \c Iterator pointing to the input to this \p transform_iterator's \c AdaptableUnaryFunction.
   *  \param f An \c AdaptableUnaryFunction used to transform the objects pointed to by \p x.
   */
  __host__ __device__
  transform_iterator(Iterator const& x, AdaptableUnaryFunction f)
      : super_t(x), m_f(f) {
  }

  /*! This explicit constructor copies the value of a given \c Iterator and creates
   *  this \p transform_iterator's \c AdaptableUnaryFunction using its null constructor.
   *
   *  \param x An \c Iterator to copy.
   */
  __host__ __device__
  explicit transform_iterator(Iterator const& x)
      : super_t(x) { }

  /*! This copy constructor creates a new \p transform_iterator from another
   *  \p transform_iterator.
   *
   *  \param other The \p transform_iterator to copy.
   */
  template<typename OtherAdaptableUnaryFunction,
           typename OtherIterator,
           typename OtherReference,
           typename OtherValue>
  __host__ __device__
  transform_iterator(const transform_iterator<OtherAdaptableUnaryFunction, OtherIterator, OtherReference, OtherValue> &other,
                            typename thrust::detail::enable_if_convertible<OtherIterator, Iterator>::type* = 0,
                            typename thrust::detail::enable_if_convertible<OtherAdaptableUnaryFunction, AdaptableUnaryFunction>::type* = 0)
      : super_t(other.base()), m_f(other.functor()) {}

  /*! Copy assignment operator copies from another \p transform_iterator.
   *  \p other The other \p transform_iterator to copy
   *  \return <tt>*this</tt>
   *
   *  \note If the type of this \p transform_iterator's functor is not copy assignable
   *        (for example, if it is a lambda) it is not an error to call this function.
   *        In this case, however, the functor will not be modified.
   *
   *        In any case, this \p transform_iterator's underlying iterator will be copy assigned.
   */
  __host__ __device__
  transform_iterator &operator=(const transform_iterator &other)
  {
    return do_assign(other,
                     // XXX gcc 4.2.1 crashes on is_copy_assignable; just assume the functor is assignable as a WAR
#if (THRUST_HOST_COMPILER == THRUST_HOST_COMPILER_GCC) && (THRUST_GCC_VERSION <= 40201)
                     thrust::detail::true_type()
#else
                     typename thrust::detail::is_copy_assignable<AdaptableUnaryFunction>::type()
#endif // THRUST_HOST_COMPILER
                     );
  }

  /*! This method returns a copy of this \p transform_iterator's \c AdaptableUnaryFunction.
   *  \return A copy of this \p transform_iterator's \c AdaptableUnaryFunction.
   */
  __host__ __device__
  AdaptableUnaryFunction functor() const
  { return m_f; }

  /*! \cond
   */
 private:
  __host__ __device__
  transform_iterator &do_assign(const transform_iterator &other, thrust::detail::true_type)
  {
    super_t::operator=(other);

    // do assign to m_f
    m_f = other.functor();

    return *this;
  }

  __host__ __device__
  transform_iterator &do_assign(const transform_iterator &other, thrust::detail::false_type)
  {
    super_t::operator=(other);

    // don't assign to m_f

    return *this;
  }

  //__thrust_hd_warning_disable__
  __host__ __device__
  typename super_t::reference dereference() const
  {
    return m_f(thrust::raw_reference_cast(*this->base()));
  }

  // tag this as mutable per Dave Abrahams in this thread:
  // http://lists.boost.org/Archives/boost/2004/05/65332.php
  mutable AdaptableUnaryFunction m_f;

  /*! \endcond
   */
}; // end transform_iterator


/*! \p make_transform_iterator creates a \p transform_iterator
 *  from an \c Iterator and \c AdaptableUnaryFunction.
 *
 *  \param it The \c Iterator pointing to the input range of the
 *            newly created \p transform_iterator.
 *  \param fun The \c AdaptableUnaryFunction used to transform the range pointed
 *             to by \p it in the newly created \p transform_iterator.
 *  \return A new \p transform_iterator which transforms the range at
 *          \p it by \p fun.
 *  \see transform_iterator
 */
template <class AdaptableUnaryFunction, class Iterator>
inline __host__ __device__
transform_iterator<AdaptableUnaryFunction, Iterator>
make_transform_iterator(Iterator it, AdaptableUnaryFunction fun)
{
  return transform_iterator<AdaptableUnaryFunction, Iterator>(it, fun);
} // end make_transform_iterator

} // end namespace v2

} // end namespace thrust



namespace fmmtl {

namespace detail {

template <typename R>
struct is_range_impl {

  template <typename T>
  static auto test(T r)
      -> decltype(*std::begin(r), void(),
                  *std::end(r), void(),
                  std::true_type());

  static auto test(...)
      -> std::false_type;

  using type = decltype(test(std::declval<R>()));
};

} // end namespace detail

template <class R>
struct is_range : detail::is_range_impl<R>::type {};

template <class R>
using iterator_type = decltype(std::begin(std::declval<R>()));


namespace range_detail {

// Poor man's quickie iterator_range
template <class Iterator>
struct iterator_range {
  iterator_range(Iterator _a, Iterator _b) : a(_a), b(_b) {}

  Iterator begin() const { return a; }
  Iterator end()   const { return b; }

  Iterator a, b;
};


template <class F, class R>
struct transformed_range
    : iterator_range<thrust::v2::transform_iterator<F, iterator_type<R> > >
{
 private:
  using iterator = thrust::v2::transform_iterator<F, iterator_type<R>>;
  using base = iterator_range<iterator>;

 public:

  transformed_range(F f, R& r)
      : base(iterator(std::begin(r), f),
             iterator(std::end(r),   f)) {
  }
};

template <class T>
struct transform_holder {
  transform_holder(const T& _t) : t(_t) {}
  T t;
};

template <class Range, class UnaryFunction>
inline transformed_range<UnaryFunction,Range>
operator|(Range& r, const transform_holder<UnaryFunction>& f) {
  return transformed_range<UnaryFunction,Range>(f.t, r);
}

template <class Range, class UnaryFunction>
inline transformed_range<UnaryFunction, const Range>
operator|(const Range& r, const transform_holder<UnaryFunction>& f) {
  return transformed_range<UnaryFunction, const Range>(f.t, r);
}

} // end namespace range_detail

template <typename F>
inline range_detail::transform_holder<F>
transformed(const F& f) {
  return {f};
}

template<class UnaryFunction, class Range>
inline range_detail::transformed_range<UnaryFunction, Range>
transform(Range& rng, UnaryFunction fn) {
  return range_detail::transformed_range<UnaryFunction, Range>(fn, rng);
}

template<class UnaryFunction, class Range>
inline range_detail::transformed_range<UnaryFunction, const Range>
transform(const Range& rng, UnaryFunction fn) {
  return range_detail::transformed_range<UnaryFunction, const Range>(fn, rng);
}

}
