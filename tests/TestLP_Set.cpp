#include "utils/LP_Set.hxx"
#include "utils/Random.h"

#include <unordered_set>

#include <gtest/gtest.h>

TEST(LP_Set, InsertContains) {

  LP_Set set(120);
  EXPECT_EQ(set.capacity(), 128);

  std::unordered_set<uint> cmp;
  auto distribution = std::uniform_int_distribution<uint>(0, std::numeric_limits<uint>::max());
  auto dice = std::bind(distribution, startGen);

  for(uint i = 0; i < 60; ++i){
    uint n = dice();
    EXPECT_EQ(set.insert(n), cmp.insert(n).second);
  }

  EXPECT_EQ(set.size(), 60);

  for(uint i : cmp) {
    EXPECT_TRUE(set.contains(i));
  }

  for(uint i : cmp) {
    EXPECT_FALSE(set.insert(i));
  }

  for(uint i = 0; i < 60; ++i){
    uint n= dice();
    EXPECT_EQ(set.contains(n), cmp.count(n));
  }
}

TEST(LP_Set, Merge) {

  LP_Set a(120);
  EXPECT_EQ(a.capacity(), 128);

  LP_Set b(120);
  EXPECT_EQ(b.capacity(), 128);

  std::unordered_set<uint> cmp;
  auto distribution = std::uniform_int_distribution<uint>(0, std::numeric_limits<uint>::max());
  auto dice = std::bind(distribution, startGen);

  for(uint i = 0; i < 60; ++i){
    uint n;
    while(n = dice(), !cmp.insert(n).second);
    EXPECT_TRUE(a.insert(n));
    while(n = dice(), !cmp.insert(n).second);
    EXPECT_TRUE(b.insert(n));
  }

  EXPECT_EQ(a.size(), 60);
  EXPECT_EQ(b.size(), 60);

  for(uint i : cmp) {
    EXPECT_TRUE(a.contains(i) | b.contains(i));
  }

  a.merge(std::move(b));

  for(uint i : cmp) {
    EXPECT_TRUE(a.contains(i));
  }
}
