#include <gtest/gtest.h>

#include "ndarray.hpp"

TEST(Generic, Generic) {
  auto *arr = GOL::NDArray<int, 3, GOL::Indexer<3,1>>::make({3,2,1});
  ASSERT_EQ( ((*arr)[{1,2,3}]), 32);
}
