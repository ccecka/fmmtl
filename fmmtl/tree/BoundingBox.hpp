#pragma once
/** @file BoundingBox.hpp
 * @brief Define the BoundingBox class for ND bounding boxes. */

#include <iostream>
#include <algorithm>
#include <cmath>

#include "fmmtl/meta/dimension.hpp"

namespace fmmtl {

/** @class BoundingBox
 * @brief Class representing ND bounding boxes.
 *
 * A BoundingBox is a ND volume. Its fundamental operations are contains(),
 * which tests whether a point is in the volume, and operator+=(), which
 * extends the volume as necessary to ensure that the volume contains a point.
 *
 * BoundingBoxes are implemented as boxes -- ND rectangular cuboids -- whose
 * sides are aligned with the principal axes.
 *
 * The point_type of a BoundingBox satisfies the concept
 * concept point_type {
 *   point_type();                             // Default constructible
 *   point_type(const point_type&);            // Copy constructible
 *   point_type& operator=(const point_type&); // Assignable
 *   value_type& operator[](unsigned i);       // mutable access to coordinate i
 * };
 * and value_type is comparable and assignable.
 */
template <typename POINT,
          unsigned DIM = fmmtl::dimension<POINT>::value>
class BoundingBox {
 public:
  static_assert(DIM >= 1, "BoundingBox point_type must have DIM >= 1");
  typedef POINT point_type;

  /** Construct an empty bounding box. */
  BoundingBox()
      : empty_(true), min_(), max_() {
  }
  /** Construct the minimal bounding box containing @a p.
   * @post contains(@a p) && min() == @a p && max() == @a p */
  explicit BoundingBox(const point_type& p)
      : empty_(false), min_(p), max_(p) {
  }
  /** Construct the minimal bounding box containing @a p1 and @a p2.
   * @post contains(@a p1) && contains(@a p2) */
  BoundingBox(const point_type& p1, const point_type& p2)
      : empty_(false), min_(p1), max_(p1) {
    *this |= p2;
  }
  /** Construct a bounding box containing the points in [first, last). */
  template <typename IT>
  BoundingBox(IT first, IT last)
      : empty_(true), min_(), max_() {
    insert(first, last);
  }

  /** Test if the bounding box is empty (contains no points). */
  bool empty() const {
    return empty_;
  }

  /** Test if the bounding box is nonempty.
   *
   * This function lets you write code such as "if (b) { ... }" or
   * "if (box1 & box2) std::cout << "box1 and box2 intersect\n". */
  operator bool() const {
    return empty();
  }

  /** Return the minimum corner of the bounding box.
   * @post empty() || contains(min())
   * @note An empty box has min() == point_type(). */
  const point_type& min() const {
    return min_;
  }

  /** Return the maximum corner of the bounding box.
   * @post empty() || contains(max())
   * @note An empty box has max() == point_type(). */
  const point_type& max() const {
    return max_;
  }

  /** Test if point @a p is in the bounding box. */
  bool contains(const point_type& p) const {
    if (empty())
      return false;
    for (unsigned i = 0; i != DIM; ++i)
      if (p[i] < min_[i] || p[i] > max_[i])
        return false;
    return true;
  }

  /** Test if @a b is entirely within this bounding box.
   * @returns true if all @a p with @a b.contains(@a p) implies contains(@a p) */
  bool contains(const BoundingBox& b) const {
    if (empty() || b.empty())
      return false;
    for (unsigned i = 0; i != DIM; ++i)
      if (b.min_[i] < min_[i] || b.min_[i] > max_[i] ||
          b.max_[i] < min_[i] || b.max_[i] > max_[i])
        return false;
    return true;
  }

  /** Test if @a b intersects this bounding box.
   * @returns true if there exists @a p such that
   *            contains(@a p) && b.contains(@a p) */
  bool intersects(const BoundingBox& b) const {
    if (empty() || b.empty())
      return false;
    for (unsigned i = 0; i != DIM; ++i)
      if (b.min_[i] > max_[i] || b.max_[i] < min_[i])
        return false;
    return true;
  }

  /** Extend the bounding box to contain @a p.
   * @post contains(@a p) is true
   * @post For all @a x with old contains(@a x),
             then new contains(@a x) is true. */
  BoundingBox& operator|=(const point_type& p) {
    if (empty()) {
      empty_ = false;
      min_ = max_ = p;
    } else {
      for (unsigned i = 0; i != DIM; ++i) {
        if (p[i] < min_[i])  min_[i] = p[i];
        if (p[i] > max_[i])  max_[i] = p[i];
      }
    }
    return *this;
  }

  /** Extend the bounding box to contain @a b.
   * @post contains(@a b) || @a b.empty()
   * @post For all @a x with old contains(@a x) or @a b.contains(@a x),
   *         then new contains(@a x) is true. */
  BoundingBox& operator|=(const BoundingBox& b) {
    if (!b.empty())
      (*this |= b.min()) |= b.max();
    return *this;
  }

  /** Extend the bounding box to contain the points in [first, last).
   * @post For all @a p in [@a first, @a last), contains(@a p) is true.
   * @post For all @a x with old contains(@a x),
   *         then new contains(@a x) is true. */
  template <typename IT>
  BoundingBox& insert(IT first, IT last) {
    for ( ; first != last; ++first)
      *this |= *first;
    return *this;
  }

  /** Intersect this bounding box with another bounding box @a b.
   * @post For all @a x with old contains(@a x) and @a b.contains(@a x),
   *         then new contains(@a x) is true. */
  BoundingBox& operator&=(const BoundingBox& b) {
    if (!intersects(b))
      return clear();
    for (unsigned i = 0; i != DIM; ++i) {
      if (min_[i] < b.min_[i])  min_[i] = b.min_[i];
      if (max_[i] > b.max_[i])  max_[i] = b.max_[i];
    }
    return *this;
  }

  /** Clear the bounding box.
   * @post empty() */
  BoundingBox& clear() {
    empty_ = true;
    min_ = max_ = point_type();
    return *this;
  }

  /** Write a BoundingBox to an output stream.
   *
   * An empty bounding box is written as "[]".
   * A nonempty bounding box is written as "[min:max] (dim)".
   */
  friend inline std::ostream& operator<<(std::ostream& s,
                                         const BoundingBox<point_type>& b) {
    if (b.empty())
      return s << '[' << ']';

    s << '[' << b.min_[0];
    for (unsigned i = 1; i != DIM; ++i)  s << ", " << b.min_[1];
    s << " : " << b.max_[0];
    for (unsigned i = 1; i != DIM; ++i)  s << ", " << b.max_[1];
    s << "] (" << b.max_[0] - b.min_[0];
    for (unsigned i = 1; i != DIM; ++i)  s << ", " << b.max_[1] - b.min_[i];
    return s << ')';
  }

 private:
  bool empty_;
  point_type min_;
  point_type max_;
};


/** Return a bounding box that contains @a b and @a p. */
template <typename P, unsigned D>
BoundingBox<P,D> operator|(BoundingBox<P,D> b,
                           const P& p) {
  return b |= p;
}
/** Return the union of @a b1 and @a b2. */
template <typename P, unsigned D>
BoundingBox<P,D> operator|(BoundingBox<P,D> b1,
                           const BoundingBox<P,D>& b2) {
  return b1 |= b2;
}
/** Return the intersection of @a b1 and @a b2. */
template <typename P, unsigned D>
BoundingBox<P,D> operator&(BoundingBox<P,D> b1,
                           const BoundingBox<P,D>& b2) {
  return b1 &= b2;
}

} // end namespace fmmtl
