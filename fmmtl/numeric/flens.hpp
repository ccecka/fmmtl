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

template <typename T>
using matrix = GeMatrix<FullStorage<T> >;

template <typename T>
using vector = DenseVector<Array<T> >;

}

//////////////////////////////
// Custom function adaptors //
//////////////////////////////

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
auto frobenius_norm(const T& t)
    -> decltype(flens::blas::asum(t)) {
  return flens::blas::asum(t);
}
