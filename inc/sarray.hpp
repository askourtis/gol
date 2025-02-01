#pragma once

#include "defs.hpp"

namespace GOL {
  template<typename T, size_t N>
  class SArray {
  public:
    T& operator[](index_t i) {
      return m_data[i];
    }

    T const& operator[](index_t i) const {
      return m_data[i];
    }
  public:
    T m_data[N];
  };
} // namespace GOL