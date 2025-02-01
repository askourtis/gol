#pragma once

#include "sarray.hpp"
#include "defs.hpp"

namespace GOL
{
  template<class T, size_t K>
  constexpr T product(SArray<T, K> const& arr) {
    T ret = 1;
    for (index_t i = 0; i < K; ++i) {
      ret *= arr[i];
    }
    return ret;
  }

} // namespace GOL
