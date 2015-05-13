#include "utils/BlockedVector.hxx"
#include <gtest/gtest.h>

TEST(BlockedVector, OneBlock) {

  BlockedVector<uint> test;

  for(uint i = 0; i < 10; ++i){
    test.emplace_back(i);
  }

  uint i = 0;
  for(auto it = test.begin(); it != test.end(); ++it){
    EXPECT_EQ(*it, i++);
  }

}

TEST(BlockedVector, OneBlockConst) {

  BlockedVector<uint> test;

  for(uint i = 0; i < 10; ++i){
    test.emplace_back(i);
  }

  uint i = 0;
  const BlockedVector<uint> & cTest = test;
  for(auto it = cTest.begin(); it != cTest.end(); ++it){
    EXPECT_EQ(*it, i++);
  }

}

TEST(BlockedVector, TwoBlocks) {

  BlockedVector<uint> a;
  BlockedVector<uint> b;

  int i;
  for(i = 0; i < 10; ++i){
    a.emplace_back(i);
  }

  for(; i < 20; ++i){
    b.emplace_back(i);
  }

  a.append(std::move(b));

  i = 0;
  for(auto it = a.begin(); it != a.end(); ++it){
    EXPECT_EQ(*it, i++);
  }
  EXPECT_EQ(i, 20);

  --i;
  for(auto it = --a.end(); it != a.begin(); --it){
    EXPECT_EQ(*it, i--);
  }
  EXPECT_EQ(i, 0);
}
