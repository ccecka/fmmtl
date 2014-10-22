#pragma once

namespace fmmtl {

// Very simple identity functor
struct identity {
  template <typename T>
  constexpr auto operator()(T&& v) const noexcept
      -> decltype(std::forward<T>(v)) {
    return std::forward<T>(v);
  }
};

}
