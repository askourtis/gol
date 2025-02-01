#include <gtest/gtest.h>

#include "ndarray.hpp"
#include "allocator.hpp"

TEST(Generic, Generic) {

  char buffer[1024];

  GOL::MonotonicPreAllocator alloc {buffer, sizeof(buffer)};

  auto &arr = GOL::make_array<int, 3>(alloc, {3,2,1});



  ASSERT_EQ((arr[GOL::DefaultIndexer<3>{1,2,3}]) , 0);
}
