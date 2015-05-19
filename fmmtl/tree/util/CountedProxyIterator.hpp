#pragma once

#include <thrust/iterator/iterator_adaptor.h>
#include <thrust/iterator/counting_iterator.h>

namespace fmmtl {

/** CountedProxyIterator
 * @brief A random-access iterator for indexed proxy objects.
 * A counted iterator over a type construted with an index-pointer pair.
 * @tparam T       The value object that supports construction: T(I, Friend*)
 * @tparam Friend  The friend class and pointer type
 * @tparam I       The index type
 *
 * Note: Since this dereferences to a value rather than a reference,
 *       it does not fully satisfy the random-access-iterator concept. Thus,
 *       this should not be implemented with thrust::transform_iterator.
 */
template <typename T, typename Friend, typename I = std::size_t>
struct CountedProxyIterator
    : thrust::iterator_adaptor<CountedProxyIterator<T,Friend,I>, // Derived
                              thrust::counting_iterator<I>,      // BaseType
                              T,                                 // Value
                              thrust::use_default,               // System(thrust)
                              std::random_access_iterator_tag,   // Traversal
                              T>                                 // Reference
{
  typedef I size_type;
  //! Construct an invalid iterator
  CountedProxyIterator() {}
  //! The index of this iterator
  size_type index() const {
    return *(this->base());
  }

  // private:
  //friend Friend;
  CountedProxyIterator(I idx, Friend* p)
      : CountedProxyIterator::iterator_adaptor(thrust::counting_iterator<I>(idx)),
        p_(p) {
  }
 private:
  Friend* p_;
  friend class thrust::iterator_core_access;
  T dereference() const {
    return T(index(), p_);
  }
};

} // end namespace fmmtl
