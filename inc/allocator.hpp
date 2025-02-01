#pragma once

#include "defs.hpp"

namespace GOL
{
  class MonotonicPreAllocator {
  public:
    MonotonicPreAllocator(void *ptr, size_t size)
      : m_ptr{ptr}, m_size{size}, m_offset{0}
    {}

    void* allocate(size_t size) {
      if (m_offset + size > m_size) {
        return nullptr;
      }
      void *ret = reinterpret_cast<char*>( m_ptr ) + m_offset;
      m_offset += size;
      return ret;
    }

    void deallocate(void* ptr, size_t size) {
      // Empty
    }

    template<class T, class ... Args>
    void construct(T *obj, Args&& ... args) {
      new (obj) T(std::forward<Args>(args)...);
    }

    template<class T>
    void destroy(T* ptr) {
      ptr->~T();
    }

  private:
    void   *m_ptr;
    size_t  m_size;
    size_t  m_offset;
  };


  class AutomaticAllocator {
  public:
    void* allocate(size_t size) {
      return nullptr;
    }

    void deallocate(void* ptr, size_t size) {
    }

    template<class T, class ... Args>
    void construct(T *obj, Args&& ... args) {
      new (obj) T(std::forward<Args>(args)...);
    }

    template<class T>
    void destroy(T* ptr) {
      ptr->~T();
    }
  };
}