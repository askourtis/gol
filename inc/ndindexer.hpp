#pragma once

#include "nddefs.hpp"
#include "meta.hpp"


namespace GOL {

  template<size_t ND, class ... IndexWrappers>
  class Indexer : public NDIndex<ND> {
    template<index_t ... Is>
    NDIndex<ND> wrapImpl(NDBoundary<ND> const& boundary, MetaIndexer<Is...>) const {
      return { IndexWrappers::wrap(this->m_data[Is], boundary[Is])... };
    }
  public:
    index_t toIndex(NDBoundary<ND> const& boundary) const;

    NDIndex<ND> wrap(NDBoundary<ND> const& boundary) const {
      return wrapImpl(boundary, make_meta_indexer<sizeof...(IndexWrappers)>());
    }
  };

  template<size_t ND, class ... IndexWrappers>
  class RowMajorIndexer : public Indexer<ND, IndexWrappers...> {
  public:
    index_t toIndex(NDBoundary<ND> const& boundary) const {
      auto idx = this->wrap(boundary);
      index_t ret = 0;
      for (size_t i = 0; i < ND; ++i) {
        ret *= boundary[i];
        ret += idx[i];
      }
      return ret;
    }
  };

  class DefaultIndexWrapper {
  public:
    static constexpr index_t wrap(index_t idx, index_t boundary) {
      return ((idx % boundary) + boundary) % boundary;
    }
  };

  template<size_t ND>
  using DefaultIndexer = RowMajorIndexer<ND, DefaultIndexWrapper, DefaultIndexWrapper, DefaultIndexWrapper>;
} // namespace GOL