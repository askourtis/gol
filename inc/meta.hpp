#pragma once

#include "defs.hpp"
namespace GOL
{

  template<index_t ... I>
  class MetaIndexer {
  public:
    template<index_t ... J>
    constexpr MetaIndexer<I..., J...> combine(MetaIndexer<J...>) {
      return {};
    }
  };

  template<size_t N>
  constexpr auto make_meta_indexer(void) {
    return make_meta_indexer<N-1>().combine(MetaIndexer<N-1>{});
  }

  template<>
  constexpr auto make_meta_indexer<0>(void) {
    return MetaIndexer<>{};
  }


  template<bool B>
  struct MetaEnableIf {
    using Type = void;
  };

  template<>
  struct MetaEnableIf<false> {
  };

  template<bool B>
  using MetaEnableIf_t = typename MetaEnableIf<B>::Type;

} // namespace GOL
