#pragma once

#include <array>
#include <utility>


namespace GOL {
  using size_t = int;
  using index_t = int;

  enum class SpaceFillingCurve : int {
    ROW_MAJOR,
    COL_MAJOR,
    Z_ORDER
  };

  template<SpaceFillingCurve T>
  class SFCSelector {};

  template<size_t CND>
  class SizeSelector {};

  template<size_t ND>
  class NDIndex {
    template<SpaceFillingCurve SFC>
    constexpr index_t toOneDimensionalIndexImpl(SFCSelector<SFC> const&);

    constexpr index_t toOneDimensionalIndexImpl(SFCSelector<SpaceFillingCurve::COL_MAJOR> const& sfcSelector) {
      return m_dimensionIndex + m_dimensionSize * m_opaqueIndex.toOneDimensionalIndex(sfcSelector);
    }

    constexpr index_t toOneDimensionalIndexImpl(SFCSelector<SpaceFillingCurve::ROW_MAJOR> const& sfcSelector) {
      return this->toOneDimensionalIndexImpl(sfcSelector, SizeSelector<1>{});
    }




    template<size_t CND>
    constexpr index_t toOneDimensionalIndexImpl(SFCSelector<SpaceFillingCurve::ROW_MAJOR> const& sfcSelector, SizeSelector<CND> const& szSelector) {
      NDIndex<CND> const *const cast_self = (NDIndex<CND> const*)this;
      return cast_self->m_dimensionIndex + cast_self->m_dimensionSize * this->toOneDimensionalIndexImpl(sfcSelector, SizeSelector<CND+1>{});
    }

    constexpr index_t toOneDimensionalIndexImpl(SFCSelector<SpaceFillingCurve::ROW_MAJOR> const& sfcSelector, SizeSelector<ND> const& szSelector) {
      return this->m_dimensionIndex;
    }

  public:
    constexpr NDIndex(std::array<std::pair<index_t, size_t>, ND> const arr) :
      m_dimensionIndex{ arr[ND-1].first }, m_dimensionSize{ arr[ND-1].second }
    {
      //Empty
    }

    template<SpaceFillingCurve SFC>
    constexpr index_t toOneDimensionalIndex(void) {
      return this->toOneDimensionalIndexImpl(SFCSelector<SFC>{});
    }

  public:
    NDIndex<ND-1> m_opaqueIndex;
    index_t       m_dimensionIndex;
    size_t        m_dimensionSize;
  };

  template<>
  class NDIndex<0> {
  public:
    template<SpaceFillingCurve SFC>
    constexpr index_t toOneDimensionalIndex(void) {
      return 0;
    }
  };




}