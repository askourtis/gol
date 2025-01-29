#pragma once



namespace GOL
{
  using size_t = int;
  using index_t = int;


  template<size_t ND>
  class NDIndexer {};


  template<class T, size_t ND, class Indexer>
  class NDArray {
    using NDIndex = std::array<index_t, ND>;
    using NDBoundary = std::array<size_t, ND>;


  public:
    static NDArray* make(NDBoundary const& boundary) {
      void *raw = std::malloc( sizeof(NDArray) + sizeof(T) * 555 );
      NDArray* cast = reinterpret_cast<NDArray*>(raw);

      cast->m_boundary = boundary;
      return reinterpret_cast<NDArray*>(raw);
    }

  public:
    constexpr T& operator[](NDIndex const& idx) {
      return m_arr[ m_indexer(idx, m_boundary) ];
    }

  private:
    Indexer    m_indexer;
    NDBoundary m_boundary;
    T          m_arr[];
  };

  template<size_t ND, int T>
  class Indexer {
  public:
    index_t operator()(std::array<index_t, ND> const& idx, std::array<index_t, ND> const& boundary) {
      return 123;
    }
  };

} // namespace GOL
