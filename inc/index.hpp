#pragma once

namespace GOL {
  using size_t = int;
  using index_t = int;

  enum class SpaceFillingCurve : int {
    ROW_MAJOR,
    COL_MAJOR,
    Z_ORDER
  };

  template<size_t ND>
  class NDIndex {
  public:
    NDIndex<ND-1> m_opaqueIndex;
    index_t m_dimensionIndex;
    size_t m_dimensionSize;
  };

  template<>
  class NDIndex<0> {};


  template<SpaceFillingCurve SFC, size_t ND>
  constexpr index_t reduceDimensions(NDIndex<ND> const& ndIndex);

  template<SpaceFillingCurve SFC>
  constexpr index_t reduceDimensions(NDIndex<0> const& ndIndex) {
    return 0;
  }

  template<size_t ND>
  constexpr index_t reduceDimensions<SpaceFillingCurve::ROW_MAJOR>(NDIndex<ND> const& ndIndex) {
    return reduceDimensions<SpaceFillingCurve::ROW_MAJOR>(ndIndex.m_opaqueIndex) + ndIndex.m_opaqueIndex.m_dimensionSize*ndIndex.m_dimensionIndex;
  }
}