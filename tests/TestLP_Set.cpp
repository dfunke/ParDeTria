#include "utils/LP_Set.hxx"
#include "utils/LP_Map.hxx"
#include "utils/LP_MultiMap.hxx"

#include <unordered_set>
#include <unordered_map>

#include <gtest/gtest.h>

TEST(Concurrent_LP_Set, InsertContains) {

    Concurrent_LP_Set set(120);
    EXPECT_EQ(set.capacity(), 128);

    std::unordered_set<uint> cmp;
    auto distribution = std::uniform_int_distribution<uint>(0, std::numeric_limits<uint>::max());
    auto dice = std::bind(distribution, startGen);

    for (uint i = 0; i < 60; ++i) {
        uint n = dice();
        EXPECT_EQ(set.insert(n), cmp.insert(n).second);
    }

    EXPECT_EQ(set.size(), 60);

    for (uint i : cmp) {
        EXPECT_TRUE(set.contains(i));
    }

    for (uint i : cmp) {
        EXPECT_FALSE(set.insert(i));
    }

    for (uint i = 0; i < 60; ++i) {
        uint n = dice();
        EXPECT_EQ(set.contains(n), cmp.count(n));
    }
}

TEST(Concurrent_LP_Set, Merge) {

    Concurrent_LP_Set a(120);
    EXPECT_EQ(a.capacity(), 128);

    Concurrent_LP_Set b(120);
    EXPECT_EQ(b.capacity(), 128);

    std::unordered_set<uint> cmp;
    auto distribution = std::uniform_int_distribution<uint>(0, std::numeric_limits<uint>::max());
    auto dice = std::bind(distribution, startGen);

    for (uint i = 0; i < 60; ++i) {
        uint n;
        while (n = dice(), !cmp.insert(n).second);
        EXPECT_TRUE(a.insert(n));
        while (n = dice(), !cmp.insert(n).second);
        EXPECT_TRUE(b.insert(n));
    }

    EXPECT_EQ(a.size(), 60);
    EXPECT_EQ(b.size(), 60);

    for (uint i : cmp) {
        EXPECT_TRUE(a.contains(i) | b.contains(i));
    }

    a.unsafe_merge(std::move(b));

    for(uint i : cmp) {
      EXPECT_TRUE(a.contains(i));
    }
}

TEST(LP_Set, InsertContains) {

    LP_Set set(120);
    EXPECT_EQ(set.capacity(), 128);

    std::unordered_set<uint> cmp;
    auto distribution = std::uniform_int_distribution<uint>(0, std::numeric_limits<uint>::max());
    auto dice = std::bind(distribution, startGen);

    for (uint i = 0; i < 60; ++i) {
        uint n = dice();
        EXPECT_EQ(set.insert(n), cmp.insert(n).second);
    }

    EXPECT_EQ(set.size(), 60);

    for (uint i : cmp) {
        EXPECT_TRUE(set.contains(i));
    }

    for (uint i : cmp) {
        EXPECT_FALSE(set.insert(i));
    }

    for (uint i = 0; i < 60; ++i) {
        uint n = dice();
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

    for (uint i = 0; i < 60; ++i) {
        uint n;
        while (n = dice(), !cmp.insert(n).second);
        EXPECT_TRUE(a.insert(n));
        while (n = dice(), !cmp.insert(n).second);
        EXPECT_TRUE(b.insert(n));
    }

    EXPECT_EQ(a.size(), 60);
    EXPECT_EQ(b.size(), 60);

    for (uint i : cmp) {
        EXPECT_TRUE(a.contains(i) | b.contains(i));
    }

    a.merge(std::move(b));

    for(uint i : cmp) {
        EXPECT_TRUE(a.contains(i));
    }
}

TEST(Concurrent_LP_Map, InsertContains) {

    Concurrent_LP_Map map(120);
    EXPECT_EQ(map.capacity(), 128);

    std::unordered_map<uint, uint> cmp;
    auto distribution = std::uniform_int_distribution<uint>(0, std::numeric_limits<uint>::max());
    auto dice = std::bind(distribution, startGen);

    for (uint i = 0; i < 60; ++i) {
        auto n = std::make_pair(dice(), dice());
        EXPECT_EQ(map.insert(n), cmp.insert(n).second);
    }

    EXPECT_EQ(map.size(), 60);

    for (auto i : cmp) {
        EXPECT_TRUE(map.contains(i.first));
        EXPECT_EQ(map.get(i.first), i.second);
    }

    for (auto i : cmp) {
        EXPECT_FALSE(map.insert(i));
    }

    for (uint i = 0; i < 60; ++i) {
        uint n = dice();
        EXPECT_EQ(map.contains(n), cmp.count(n));
    }
}

TEST(Concurrent_LP_Map, Merge) {

    Concurrent_LP_Map a(120);
    EXPECT_EQ(a.capacity(), 128);

    Concurrent_LP_Map b(120);
    EXPECT_EQ(b.capacity(), 128);

    std::unordered_map<uint, uint> cmp;
    auto distribution = std::uniform_int_distribution<uint>(0, std::numeric_limits<uint>::max());
    auto dice = std::bind(distribution, startGen);

    for (uint i = 0; i < 60; ++i) {
        std::pair<uint, uint> n;
        while (n = std::make_pair(dice(), dice()), !cmp.insert(n).second);
        EXPECT_TRUE(a.insert(n));
        while (n = std::make_pair(dice(), dice()), !cmp.insert(n).second);
        EXPECT_TRUE(b.insert(n));
    }

    EXPECT_EQ(a.size(), 60);
    EXPECT_EQ(b.size(), 60);

    for (auto i : cmp) {
        EXPECT_TRUE(a.contains(i.first) | b.contains(i.first));
        if(a.contains(i.first)) {
            EXPECT_EQ(a.get(i.first), i.second);
            EXPECT_EQ(b.get(i.first), 0);
        }
        if(b.contains(i.first)) {
            EXPECT_EQ(a.get(i.first), 0);
            EXPECT_EQ(b.get(i.first), i.second);
        }
    }

    a.unsafe_merge(std::move(b));

    for(auto i : cmp) {
        EXPECT_TRUE(a.contains(i.first));
        EXPECT_EQ(a.get(i.first), i.second);
    }
}

TEST(Concurrent_LP_MultiMap, InsertContains) {

    Concurrent_LP_MultiMap map(600);
    EXPECT_EQ(map.capacity(), 1024);

    std::unordered_multimap<uint, uint> cmp;

    auto valueDistribution = std::uniform_int_distribution<uint>(0, std::numeric_limits<uint>::max());
    auto valueDice = std::bind(valueDistribution, startGen);

    auto numberDistribution = std::uniform_int_distribution<uint>(1, 10);
    auto numberDice = std::bind(numberDistribution, startGen);

    for (uint i = 0; i < 60; ++i) {
        auto key = valueDice();
        auto n = numberDice();
        for(uint k = 0; k < n; ++k) {
            auto x = std::make_pair(key, valueDice());

            EXPECT_TRUE(map.insert(x));
            cmp.insert(x);
        }
    }

    EXPECT_EQ(map.size(), cmp.size());

    for (auto i : cmp) {
        EXPECT_TRUE(map.contains(i.first));

        auto cmpRange = cmp.equal_range(i.first);
        auto mapRange = map.get(i.first);

        EXPECT_EQ(std::distance(mapRange.first, mapRange.second), std::distance(cmpRange.first, cmpRange.second));

        for(auto it = mapRange.first; it != mapRange.second; ++it){
            EXPECT_TRUE(std::find(cmpRange.first, cmpRange.second, *it) != cmpRange.second);
        }

    }
}

TEST(Concurrent_LP_MultiMap, Merge) {

    Concurrent_LP_MultiMap a(600);
    EXPECT_EQ(a.capacity(), 1024);

    Concurrent_LP_MultiMap b(600);
    EXPECT_EQ(b.capacity(), 1024);

    std::unordered_multimap<uint, uint> cmp;

    auto valueDistribution = std::uniform_int_distribution<uint>(0, std::numeric_limits<uint>::max());
    auto valueDice = std::bind(valueDistribution, startGen);

    auto numberDistribution = std::uniform_int_distribution<uint>(1, 10);
    auto numberDice = std::bind(numberDistribution, startGen);

    for (uint i = 0; i < 60; ++i) {
        auto key = valueDice();

        auto n = numberDice();
        for(uint k = 0; k < n; ++k) {
            auto x = std::make_pair(key, valueDice());

            EXPECT_TRUE(a.insert(x));
            cmp.insert(x);
        }

        n = numberDice();
        for(uint k = 0; k < n; ++k) {
            auto x = std::make_pair(key, valueDice());

            EXPECT_TRUE(b.insert(x));
            cmp.insert(x);
        }
    }

    EXPECT_EQ(a.size() + b.size(), cmp.size());

    for (auto it = cmp.begin(); it != cmp.end();) {
        auto i = *it;

        EXPECT_TRUE(a.contains(i.first) && b.contains(i.first));

        auto cmpRange = cmp.equal_range(i.first);
        auto aRange = a.get(i.first);
        auto bRange = b.get(i.first);

        EXPECT_EQ(std::distance(aRange.first, aRange.second) + std::distance(bRange.first, bRange.second), std::distance(cmpRange.first, cmpRange.second));

        for(auto it = aRange.first; it != aRange.second; ++it){
            EXPECT_TRUE(std::find(cmpRange.first, cmpRange.second, *it) != cmpRange.second);
        }

        for(auto it = bRange.first; it != bRange.second; ++it){
            EXPECT_TRUE(std::find(cmpRange.first, cmpRange.second, *it) != cmpRange.second);
        }

        it = cmpRange.second;
    }

    a.unsafe_merge(std::move(b));

    for (auto it = cmp.begin(); it != cmp.end();) {
        auto i = *it;

        EXPECT_TRUE(a.contains(i.first));

        auto cmpRange = cmp.equal_range(i.first);
        auto mapRange = a.get(i.first);

        EXPECT_EQ(std::distance(mapRange.first, mapRange.second), std::distance(cmpRange.first, cmpRange.second));

        for(auto it = mapRange.first; it != mapRange.second; ++it){
            EXPECT_TRUE(std::find(cmpRange.first, cmpRange.second, *it) != cmpRange.second);
        }

        it = cmpRange.second;
    }
}