#pragma once

#ifndef USE_CXXLAPACK
#  define USE_CXXLAPACK
#endif

#include <flens/flens.cxx>

/////////////////////
// Custom adaptors //
/////////////////////

template <typename T>
auto num_rows(T&& t)
    -> decltype(t.numRows()) {
  return t.numRows();
}

template <typename T>
auto num_cols(T&& t)
    -> decltype(t.numCols()) {
  return t.numCols();
}
