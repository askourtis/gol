#pragma once

#include <initializer_list>

#include "nddefs.hpp"
#include "ndindexer.hpp"
#include "allocator.hpp"
#include "algorithm.hpp"

namespace GOL
{
  template<class T, size_t ND, class Allocator = AutomaticAllocator>
  class NDArray {
  public:
    NDArray(Allocator const& allocator, NDBoundary<ND> const& boundary, T* arr)
      : m_allocator{allocator}, m_boundary{boundary}, m_arr{arr}
    {}

    ~NDArray(void) {
      for (size_t i = 0; i < this->size() ; ++i) {
        m_allocator.destroy(&m_arr[i]);
      }
      m_allocator.deallocate(this, this->getByteCount());
    }

    template<class Indexer>
    constexpr T& operator[](Indexer const& idx) {
      return m_arr[ idx.toIndex(m_boundary) ];
    }

    constexpr size_t size() {
      return product(m_boundary);
    }

    constexpr NDBoundary<ND> const& boundary() {
      return m_boundary;
    }

  private:
    constexpr size_t getByteCount() {
      size_t arr_sz = sizeof(T);
      for (auto const& n : m_boundary) {
        arr_sz *= n;
      }
      return arr_sz + sizeof(*this);
    }

  private:
    Allocator       m_allocator;
    NDBoundary<ND>  m_boundary;
    T              *m_arr;
  };

  template<class T, size_t ND, class Allocator>
  NDArray<T, ND, Allocator> &make_array(Allocator& alloc, NDBoundary<ND> const& boundary) {
    void *raw = alloc.allocate( sizeof(NDArray<T, ND, Allocator>) + sizeof(T) * product(boundary) );
    auto *cast = reinterpret_cast<NDArray<T, ND, Allocator>*>(raw);
    alloc.construct(cast, alloc, boundary, reinterpret_cast<T*>( cast+1 ));
    return *cast;
  }
} // namespace GOL
