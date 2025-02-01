#pragma once

#include "nddefs.hpp"


namespace GOL {

  template<size_t ND>
  class Indexer : public SArray<index_t, ND> {
  public:
    index_t toIndex(NDBoundary<ND> const& boundary) const;
  };

  template<size_t ND>
  class ColMajorIndexer : public Indexer<ND> {
  public:
    index_t toIndex(NDBoundary<ND> const& boundary) const {
      index_t ret = 0;
      for (size_t i = 0; i < ND; ++i) {
        ret *= boundary[i];
        ret += (*this)[i];
      }
      return ret;
    }
  };

  template<size_t ND, class IndexWrapper, class ... IndexWrapperArgs>
  class RowMajorIndexer : public Indexer<ND> {
  public:
    index_t toIndex(NDBoundary<ND> const& boundary) const {
      return IndexWrapper::wrap( (*this)[ND-1], boundary[ND-1] ) * boundary[ND-1] + RowMajorIndexer<ND-1, IndexWrapperArgs...>(this->data).toIndex(boundary);
    }
  };

  template<size_t ND>
  using DefaultIndexer = ColMajorIndexer<ND>;
} // namespace GOL