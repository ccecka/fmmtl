#pragma once

#ifndef ALWAYS_USE_CXXLAPACK
#  define ALWAYS_USE_CXXLAPACK
#endif

#include <flens/flens.cxx>

//////////////////////////
// Custom type adaptors //
//////////////////////////

// TODO: Full MTL4-like type options/construction

namespace flens {

/////////////////////////
// Convenience Typedef //
/////////////////////////

template <typename T>
using matrix = GeMatrix<FullStorage<T> >;

template <typename T>
using vector = DenseVector<Array<T> >;

/////////////////////////
// Convenience Methods //
/////////////////////////

template <typename T>
auto num_rows(const T& t)
    -> decltype(t.numRows()) {
  return t.numRows();
}

template <typename T>
auto num_cols(const T& t)
    -> decltype(t.numCols()) {
  return t.numCols();
}

template <typename T>
auto norm_f(const T& t)
    -> decltype(flens::blas::asum(t)) {
  return flens::blas::asum(t);
}

template <typename IndexType, typename OtherIndexType>
Range<IndexType>
operator-(const Range<IndexType>& r, OtherIndexType i) {
  return Range<IndexType>(r.firstIndex() - i, r.stride(), r.lastIndex() - i);
};

template <typename IndexType, typename OtherIndexType>
Range<IndexType>
operator+(const Range<IndexType>& r, OtherIndexType i) {
  return Range<IndexType>(r.firstIndex() + i, r.stride(), r.lastIndex() + i);
};

}
