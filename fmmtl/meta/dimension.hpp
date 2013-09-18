#pragma once

namespace fmmtl {

template <typename T>
struct dimension;

template <>
struct dimension<double> {
  const static unsigned value = 1;
};

template <>
struct dimension<float> {
  const static unsigned value = 1;
};

} // end namepsace fmmtl
