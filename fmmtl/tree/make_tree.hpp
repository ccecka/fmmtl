#pragma once

#include "Octree.hpp"

template <class Tree>
struct TreeBuilder {
  template <class PointIter,
            class Options>
  static inline Tree make(PointIter first, PointIter last, Options& opts);
};


#include "Octree.hpp"

template <class P>
template <class PointIter,
          class Options>
static inline
Octree<P>
TreeBuilder<Octree<P>>::make(PointIter first, PointIter last, Options& opts) {
  return Octree<P>(first, last, opts.ncrit());
}
