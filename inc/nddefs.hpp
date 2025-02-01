#pragma once

#include "defs.hpp"
#include "sarray.hpp"

namespace GOL
{
  template<size_t ND>
  using NDIndex = SArray<index_t, ND>;

  template<size_t ND>
  using NDBoundary = SArray<size_t, ND>;
} // namespace GOL
