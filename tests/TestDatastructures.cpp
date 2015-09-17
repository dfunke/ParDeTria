#include "datastructures/LP_Set.hxx"
#include "datastructures/LP_Map.hxx"
#include "datastructures/LP_MultiMap.hxx"
#include "datastructures/BlockedArray.hxx"
#include "datastructures/Growing_LP.hxx"
#include "datastructures/Bit_Set.hxx"
#include "datastructures/BlockedArray2.hxx"

#include <unordered_set>
#include <unordered_map>
#include <algorithm>

#include <gtest/gtest.h>

TEST(Concurrent_LP_Set, InsertContains) {

    Concurrent_LP_Set set(120);
    EXPECT_EQ(set.capacity(), 128);

    std::unordered_set<uint> cmp;
    auto distribution = std::uniform_int_distribution<uint>(1, std::numeric_limits<uint>::max());
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
    auto distribution = std::uniform_int_distribution<uint>(1, std::numeric_limits<uint>::max());
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

    a.unsafe_merge(std::move(b), std::set<uint>());

    EXPECT_EQ(a.size(), 120);

    for (uint i : cmp) {
        EXPECT_TRUE(a.contains(i));
    }
}

TEST(Concurrent_LP_Set, ConcurrentGrow) {

    GrowingHashTable<Concurrent_LP_Set> set(120);
    GrowingHashTableHandle<Concurrent_LP_Set> handle(set);
    EXPECT_EQ(handle.capacity(), 128);

    tbb::concurrent_unordered_set<uint> cmp;

    //auto distribution = std::uniform_int_distribution<uint>(1, std::numeric_limits<uint>::max());
    //auto dice = std::bind(distribution, startGen);

    tbb::parallel_for(tbb::blocked_range<uint>(1, 1e6), [&](const auto &r) {
        GrowingHashTableHandle<Concurrent_LP_Set> handle(set);
        for (auto i = r.begin(); i != r.end(); ++i) {
            handle.insert(i);
            cmp.insert(i);
        }
    });

    EXPECT_EQ(cmp.size(), handle.size());

    for (uint i : cmp) {
        EXPECT_TRUE(handle.contains(i));
    }
}

TEST(LP_Set, InsertContains) {

    LP_Set set(120);
    EXPECT_EQ(set.capacity(), 128);

    std::unordered_set<uint> cmp;
    auto distribution = std::uniform_int_distribution<uint>(1, std::numeric_limits<uint>::max());
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
    auto distribution = std::uniform_int_distribution<uint>(1, std::numeric_limits<uint>::max());
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

    for (uint i : cmp) {
        EXPECT_TRUE(a.contains(i));
    }
}

TEST(Concurrent_LP_Map, InsertContains) {

    Concurrent_LP_Map map(120);
    EXPECT_EQ(map.capacity(), 128);

    std::unordered_map<uint, uint> cmp;
    auto distribution = std::uniform_int_distribution<uint>(1, std::numeric_limits<uint>::max());
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
    auto distribution = std::uniform_int_distribution<uint>(1, std::numeric_limits<uint>::max());
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
        if (a.contains(i.first)) {
            EXPECT_EQ(a.get(i.first), i.second);
            EXPECT_EQ(b.get(i.first), 0);
        }
        if (b.contains(i.first)) {
            EXPECT_EQ(a.get(i.first), 0);
            EXPECT_EQ(b.get(i.first), i.second);
        }
    }

    a.unsafe_merge(std::move(b));

    for (auto i : cmp) {
        EXPECT_TRUE(a.contains(i.first));
        EXPECT_EQ(a.get(i.first), i.second);
    }
}

TEST(Concurrent_LP_MultiMap, InsertContains) {

    Concurrent_LP_MultiMap map(600);
    EXPECT_EQ(map.capacity(), 1024);

    std::unordered_multimap<uint, uint> cmp;

    auto valueDistribution = std::uniform_int_distribution<uint>(1, std::numeric_limits<uint>::max());
    auto valueDice = std::bind(valueDistribution, startGen);

    auto numberDistribution = std::uniform_int_distribution<uint>(1, 10);
    auto numberDice = std::bind(numberDistribution, startGen);

    for (uint i = 0; i < 60; ++i) {
        auto key = valueDice();
        auto n = numberDice();
        for (uint k = 0; k < n; ++k) {
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

        for (auto it = mapRange.first; it != mapRange.second; ++it) {
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

    auto valueDistribution = std::uniform_int_distribution<uint>(1, std::numeric_limits<uint>::max());
    auto valueDice = std::bind(valueDistribution, startGen);

    auto numberDistribution = std::uniform_int_distribution<uint>(1, 10);
    auto numberDice = std::bind(numberDistribution, startGen);

    for (uint i = 0; i < 60; ++i) {
        auto key = valueDice();

        auto n = numberDice();
        for (uint k = 0; k < n; ++k) {
            auto x = std::make_pair(key, valueDice());

            EXPECT_TRUE(a.insert(x));
            cmp.insert(x);
        }

        n = numberDice();
        for (uint k = 0; k < n; ++k) {
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

        EXPECT_EQ(std::distance(aRange.first, aRange.second) + std::distance(bRange.first, bRange.second),
                  std::distance(cmpRange.first, cmpRange.second));

        for (auto it = aRange.first; it != aRange.second; ++it) {
            EXPECT_TRUE(std::find(cmpRange.first, cmpRange.second, *it) != cmpRange.second);
        }

        for (auto it = bRange.first; it != bRange.second; ++it) {
            EXPECT_TRUE(std::find(cmpRange.first, cmpRange.second, *it) != cmpRange.second);
        }

        it = cmpRange.second;
    }

    a.unsafe_merge(std::move(b), std::set<uint>());

    for (auto it = cmp.begin(); it != cmp.end();) {
        auto i = *it;

        EXPECT_TRUE(a.contains(i.first));

        auto cmpRange = cmp.equal_range(i.first);
        auto mapRange = a.get(i.first);

        EXPECT_EQ(std::distance(mapRange.first, mapRange.second), std::distance(cmpRange.first, cmpRange.second));

        for (auto it = mapRange.first; it != mapRange.second; ++it) {
            EXPECT_TRUE(std::find(cmpRange.first, cmpRange.second, *it) != cmpRange.second);
        }

        it = cmpRange.second;
    }
}

TEST(BlockedArray, OneBlock) {

    const uint BS = 512;

    BlockedArray<uint, BS> ba; // allocates one block;

    for (uint i = 0; i < BS; ++i)
        ba[i] = i;

    for (uint i = 0; i < BS; ++i)
            EXPECT_EQ(i, ba[i]);

}

TEST(BlockedArray, TwoBlocks) {

    const uint BS = 512;

    BlockedArray<uint, BS> ba(2 * BS); // allocates two block;

    for (uint i = 0; i < 2 * BS; ++i)
        ba[i] = i;

    for (uint i = 0; i < 2 * BS; ++i)
            EXPECT_EQ(i, ba[i]);

}

TEST(BlockedArray, ReserveBlock) {

    const uint BS = 512;

    BlockedArray<uint, BS> ba; // allocates one block;

    for (uint i = 0; i < 2 * BS; ++i) {
        if (i == BS / 2)
            ba.reserve(2 * BS);

        ba[i] = i;
    }

    for (uint i = 0; i < 2 * BS; ++i)
            EXPECT_EQ(i, ba[i]);

}

TEST(Concurrent_BlockedArray, OneBlock) {

    const uint BS = 512;

    Concurrent_BlockedArray<uint, BS> ba; // allocates one block;

    for (uint i = 0; i < BS; ++i)
        ba[i] = i;

    for (uint i = 0; i < BS; ++i)
            EXPECT_EQ(i, ba[i]);

}

TEST(Concurrent_BlockedArray, TwoBlocks) {

    const uint BS = 512;

    Concurrent_BlockedArray<uint, BS> ba(2 * BS); // allocates two block;

    for (uint i = 0; i < 2 * BS; ++i)
        ba[i] = i;

    for (uint i = 0; i < 2 * BS; ++i)
            EXPECT_EQ(i, ba[i]);

}

TEST(Concurrent_BlockedArray, ReserveBlock) {

    const uint BS = 512;

    Concurrent_BlockedArray<uint, BS> ba; // allocates one block;

    for (uint i = 0; i < 2 * BS; ++i) {
        if (i == BS / 2)
            ba.reserve(2 * BS);

        ba[i] = i;
    }

    for (uint i = 0; i < 2 * BS; ++i)
            EXPECT_EQ(i, ba[i]);

}

TEST(Bit_Set, InsertContains) {

    const uint cap = 1e5;
    Bit_Set set(0, cap);
    EXPECT_GE(set.capacity(), cap);

    std::unordered_set<uint> cmp;
    auto distribution = std::uniform_int_distribution<uint>(1, cap);
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

TEST(Bit_Set, BatchInsert) {

    const uint cap = 320;

    auto setAndTest = [cap](const uint from, const uint to) {
        Bit_Set set(0, cap);
        EXPECT_GE(set.capacity(), cap);

        set.batchSet(from, to);

        std::size_t i;
        for (i = 0; i < from; ++i)
            EXPECT_TRUE(set.zero(i));
        for (; i < to; ++i)
            EXPECT_TRUE(set.one(i));
        for (; i < set.capacity(); ++i)
            EXPECT_TRUE(set.zero(i));
    };

    //within one block
    setAndTest(7, 29);

    //one entire block
    setAndTest(0, 32);
    setAndTest(287, 320);

    //two entire blocks
    setAndTest(0, 64);
    setAndTest(255, 320);

    //two broken blocks
    setAndTest(9, 53);
    setAndTest(275, 309);

    //several entire blocks
    setAndTest(0, 320);
    setAndTest(95, 320);

    //several broken blocks
    setAndTest(5, 309);
    setAndTest(105, 289);
}

TEST(Bit_Set, Merge) {


    const uint capA = 1e5;
    Bit_Set a(0, capA);
    EXPECT_GE(a.capacity(), capA);

    const uint capB = 1e6;
    Bit_Set b(0, capB);
    EXPECT_GE(b.capacity(), capB);

    std::unordered_set<uint> cmp;
    auto distribution = std::uniform_int_distribution<uint>(1, std::min(capA, capB));
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
    EXPECT_GE(a.capacity(), std::max(capA, capB));

    for (uint i : cmp) {
        EXPECT_TRUE(a.contains(i));
    }
}

TEST(Bit_Set, Resize) {


    const uint cap = 1e5;
    Bit_Set set(0, cap);
    EXPECT_GE(set.capacity(), cap);

    std::unordered_set<uint> cmp;
    auto distribution = std::uniform_int_distribution<uint>(1, cap);
    auto dice = std::bind(distribution, startGen);

    for (uint i = 0; i < 60; ++i) {
        uint n;
        while (n = dice(), !cmp.insert(n).second);
        EXPECT_TRUE(set.insert(n));
    }

    EXPECT_EQ(set.size(), 60);

    for (uint i : cmp) {
        EXPECT_TRUE(set.contains(i));
    }

    set.resize(2 * cap);
    EXPECT_GE(set.capacity(), 2 * cap);

    for (uint i = 0; i < set.capacity(); ++i) {
        EXPECT_EQ(set.contains(i), cmp.count(i));
    }
}

TEST(Bit_Set, Filter) {


    const uint cap = 1e5;
    Bit_Set set(0, cap);
    EXPECT_GE(set.capacity(), cap);

    set.batchSet(0, cap);
    EXPECT_EQ(set.ones(), cap);

    Bit_Set filter(0, cap);
    for (uint i = 1; i < cap; i += 2)
        filter.insert(i);

    set.filter(filter);

    EXPECT_EQ(set.ones(), cap / 2);
    for (uint i = 0; i < cap; ++i)
            EXPECT_EQ(set.one(i), i % 2 == 0);
}

TEST(Bit_Set, MergeFilter) {


    const uint capA = 1e5;
    Bit_Set a(0, capA);
    EXPECT_GE(a.capacity(), capA);

    a.batchSet(0, capA);
    EXPECT_EQ(a.ones(), capA);

    const uint capB = 1e6;
    Bit_Set b(0, capB);
    EXPECT_GE(b.capacity(), capB);

    b.batchSet(0, capB);
    EXPECT_EQ(b.ones(), capB);

    Bit_Set filter(0, std::max(capA, capB));
    for (uint i = 1; i < std::max(capA, capB); i += 2)
        filter.insert(i);

    a.mergeFilter(std::move(b), filter);
    EXPECT_GE(a.capacity(), std::max(capA, capB));

    EXPECT_EQ(a.ones(), std::max(capA, capB) / 2);
    for (uint i = 0; i < std::max(capA, capB); ++i)
            EXPECT_EQ(a.one(i), i % 2 == 0);
}

TEST(Concurrent_Bit_Set, InsertContains) {

    const uint cap = 1e6;
    Concurrent_Bit_Set set(0, cap);
    EXPECT_GE(set.capacity(), cap);

    tbb::concurrent_unordered_set<uint> cmp;
    auto distribution = std::uniform_int_distribution<uint>(1, cap);
    auto dice = std::bind(distribution, startGen);

    tbb::parallel_for(tbb::blocked_range<uint>(1, cap / 10), [&](__attribute__((unused)) const auto &r) {
        uint n = dice();
        set.insert(n);
        cmp.insert(n);
    });

    EXPECT_EQ(set.size(), cmp.size());

    for (uint i : cmp) {
        EXPECT_TRUE(set.contains(i));
    }

    for (uint i : cmp) {
        EXPECT_FALSE(set.insert(i));
    }

    tbb::parallel_for(tbb::blocked_range<uint>(1, cap / 10), [&](__attribute__((unused)) const auto &r) {
        uint n = dice();
        bool contained = cmp.count(n);

        EXPECT_EQ(contained, set.contains(n));
    });
}

TEST(BlockedArray2, InsertMerge) {

    const uint N = 120;
    BlockedArray2<uint> a(0, N);

    // test iterator and set array to consecutive numbers
    uint count = 0;
    for (auto &i : a) {
        i = count++;
    }
    EXPECT_EQ(N, count);

    // test for correct values
    for (uint i = 0; i < N; ++i) {
        EXPECT_EQ(i, a.unsafe_at(i));
    }

    EXPECT_THROW(a.unsafe_at(N), std::out_of_range);

    BlockedArray2<uint> b(N, 2 * N);

    // test iterator and set array to consecutive numbers
    count = 0;
    for (auto &i : b) {
        i = count++;
    }
    EXPECT_EQ(N, count);

    // test for correct values
    for (uint i = N; i < 2 * N; ++i) {
        EXPECT_EQ(i - N, b.unsafe_at(i));
    }

    EXPECT_THROW(b.unsafe_at(2 * N), std::out_of_range);
    EXPECT_THROW(b.unsafe_at(N - 1), std::out_of_range);

    a.merge(std::move(b));

    // test for correct values
    for (uint i = 0; i < 2 * N; ++i) {
        EXPECT_EQ(i % N, a.unsafe_at(i));
    }

    // test iterator
    for (const auto &i : a) {
        EXPECT_GE(i, 0);
        EXPECT_LT(i, N);
    }

    EXPECT_THROW(a.unsafe_at(2 * N), std::out_of_range);
}

TEST(BlockedArray2, ParallelRange) {

    const uint N = 1e5;
    BlockedArray2<uint> a(0, N);
    // test iterator and set array to consecutive numbers
    uint count = 0;
    for (auto &i : a) {
        i = count++;
    }
    EXPECT_EQ(N, count);

    tbb::concurrent_unordered_set<uint> seen;
    tbb::parallel_for(a.range(), [N, &seen](const auto &r) {
        for (const auto i : r) {
            EXPECT_GE(i, 0);
            EXPECT_LT(i, N);
            seen.insert(i);
        }
    });
    EXPECT_EQ(N, seen.size());

    BlockedArray2<uint> b(2 * N, 3 * N);
    // test iterator and set array to consecutive numbers
    count = 2 * N;
    for (auto &i : b) {
        i = count++;
    }
    EXPECT_EQ(3 * N, count);

    seen.clear();
    tbb::parallel_for(b.range(), [N, &seen](const auto &r) {
        for (const auto i : r) {
            EXPECT_GE(i, 2 * N);
            EXPECT_LT(i, 3 * N);
            seen.insert(i);
        }
    });
    EXPECT_EQ(N, seen.size());

    a.merge(std::move(b));

    seen.clear();
    tbb::parallel_for(a.range(), [N, &seen](const auto &r) {
        for (const auto i : r) {
            EXPECT_GE(i, 0);
            EXPECT_LT(i, 3 * N);
            seen.insert(i);
        }
    });
    EXPECT_EQ(2 * N, seen.size());
}