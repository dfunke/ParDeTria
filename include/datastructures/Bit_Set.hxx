#pragma once

#include "datastructures/LP_Set.hxx"

#include "utils/Misc.h"
#include "utils/ASSERT.h"

class Concurrent_Bit_Set;

class Bit_Set {

    friend class Concurrent_Bit_Set;

public:
    typedef uint tKeyType;

public:

    Bit_Set(std::size_t size)
            : m_array(nullptr),
              m_ones(0) {

        // Initialize cells
        m_arraySize = (size + cENTRIES - 1) / cENTRIES; //get the ceiling
        m_capacity = m_arraySize * cENTRIES;
        try {
            m_array.reset(new tKeyType[m_arraySize]()); //value init
        } catch (std::bad_alloc &e) {
            std::cerr << e.what() << std::endl;
            raise(SIGINT);
        }
    }

    Bit_Set(Bit_Set &&other)
            : m_arraySize(other.m_arraySize),
              m_capacity(other.m_capacity),
              m_array(std::move(other.m_array)),
              m_ones(other.m_ones) { }

    //Conversion
    Bit_Set(Concurrent_Bit_Set &&other);

    Bit_Set &operator=(Bit_Set &&other) {
        m_arraySize = other.m_arraySize;
        m_capacity = other.m_capacity;
        m_array = std::move(other.m_array);
        m_ones = other.m_ones;

        return *this;
    }

    bool testAndSet(const tKeyType &idx) {
        RAISE(idx / cENTRIES < m_arraySize);
        RAISE(idx < m_capacity);

        bool test = isSet(idx);
        m_array[idx / cENTRIES] |= 1 << (idx % cENTRIES);
        m_ones += !test;

        return test;
    }

    void batchSet(const tKeyType &from, const tKeyType &to) {

        /*auto print = [this] (const std::string & msg) {
            std::cout << msg << ": ";
            for(std::size_t i = 0; i < capacity(); ++i)
                std::cout << isSet(i);
            std::cout << std::endl;
        };*/

        std::size_t fromBlock = from / cENTRIES;
        std::size_t toBlock = (to - 1) / cENTRIES; // to not included

        //print("before");

        // special case: first block
        m_array[fromBlock] = (~(tKeyType) 0) << (from % cENTRIES);

        //print("first ");

        // middle block: all ones
        for (std::size_t i = fromBlock + 1; i < toBlock; ++i)
            m_array[i] = (~(tKeyType) 0);

        //print("middle");

        // special case: last block
        m_array[toBlock] = ((~(tKeyType) 0) >> (cENTRIES - (to % cENTRIES))) &
                           (toBlock == fromBlock ? m_array[toBlock] : (~(tKeyType) 0));

        m_ones = to - from;


        //print("after ");
    }

    bool isSet(const tKeyType &idx) const {

        if (idx < m_capacity) {
            RAISE(idx / cENTRIES < m_arraySize);

            return m_array[idx / cENTRIES] && (m_array[idx / cENTRIES] & (1 << (idx % cENTRIES)));
        } else {
            return false;
        }
    }

    bool zero() const { return m_ones == 0; }

    bool one() const { return m_ones == m_capacity; }

    bool zero(const std::size_t idx) const {
        return !isSet(idx);
    }

    bool one(const std::size_t idx) const {
        return isSet(idx);
    }

    std::size_t capacity() const { return m_capacity; }

    std::size_t zeros() const { return m_capacity - m_ones; }

    std::size_t ones() const { return m_ones; }

    template<class BitSet>
    void copy(const BitSet &other) {

        m_arraySize = other.m_arraySize;
        m_capacity = other.m_capacity;
        m_ones = other.m_ones;

        try {
            m_array.reset(new tKeyType[m_arraySize]); //random init
        } catch (std::bad_alloc &e) {
            std::cerr << e.what() << std::endl;
            raise(SIGINT);
        }

        for (std::size_t i = 0; i < m_arraySize; ++i) {
            m_array[i] = other.m_array[i];
        }
    }

    void resize(const std::size_t &newSize) {

        tArrayPtr oldArray = std::move(m_array);
        std::size_t oldSize = m_arraySize;

        m_arraySize = (newSize + cENTRIES - 1) / cENTRIES; //get the ceiling
        m_capacity = m_arraySize * cENTRIES;

        try {
            m_array.reset(new tKeyType[m_arraySize]); //random init
        } catch (std::bad_alloc &e) {
            std::cerr << e.what() << std::endl;
            raise(SIGINT);
        }

        std::size_t i;
        for (i = 0; i < oldSize; ++i) {
            m_array[i] = oldArray[i];
        }
        for (; i < m_arraySize; ++i) {
            m_array[i] = (tKeyType) 0;
        }
    }

    template<class BitSet>
    void merge(BitSet &&other) {

        std::size_t bSize;
        tArrayPtr bArray;
        if (m_capacity >= other.m_capacity) {
            bSize = other.m_arraySize;
            bArray = std::move(other.m_array);
        } else {
            bSize = m_arraySize;
            bArray = std::move(m_array);

            m_arraySize = other.m_arraySize;
            m_capacity = other.m_capacity;
            m_array = std::move(other.m_array);
        }

        m_ones = 0;
        for (std::size_t i = 0; i < bSize; ++i) {
            m_array[i] |= bArray[i];
            m_ones += popcount(m_array[i]);
        }

    }

    template<class BitSet>
    void filter(const BitSet &filter) {

        m_ones = 0;
        for (std::size_t i = 0; i < std::min(m_arraySize, filter.m_arraySize); ++i) {
            m_array[i] &= ~filter.m_array[i];
            m_ones += popcount(m_array[i]);
        }
    }

    template<class BitSet, class Filter>
    void mergeFilter(BitSet &&other, const Filter &filter) {

        std::size_t fSize = filter.m_arraySize;

        std::size_t bSize;
        tArrayPtr bArray;
        if (m_capacity >= other.m_capacity) {
            bSize = other.m_arraySize;
            bArray = std::move(other.m_array);
        } else {
            bSize = m_arraySize;
            bArray = std::move(m_array);

            m_arraySize = other.m_arraySize;
            m_capacity = other.m_capacity;
            m_array = std::move(other.m_array);
        }

        m_ones = 0;
        std::size_t i;
        for (i = 0; i < bSize; ++i) {
            m_array[i] = (m_array[i] | bArray[i]) & ~(i < fSize ? filter.m_array[i] : (tKeyType) 0);
            m_ones += popcount(m_array[i]);
        }
        for (; i < fSize; ++i) {
            m_array[i] &= ~filter.m_array[i];
            m_ones += popcount(m_array[i]);
        }

    }

    // compat definitions
    bool insert(const tKeyType &key) { return !testAndSet(key); }

    bool contains(const tKeyType &key) const { return isSet(key); }

    std::size_t count(const tKeyType &key) const { return contains(key); }

    bool empty() const { return zero(); };

    bool empty(const std::size_t idx) const { return zero(idx); }

    tKeyType at(const std::size_t idx) const { return idx * isSet(idx); }

    std::size_t size() const { return ones(); }


public:
    typedef _detail::iterator<Bit_Set, tKeyType> iterator;
    typedef _detail::range_type<Bit_Set, iterator> range_type;

    iterator begin() const {
        return iterator(*this, 0, true);
    }

    iterator end() const {
        return iterator(*this, m_capacity, false);
    }

    range_type range() const {
        return range_type(*this);
    }

private:
    typedef std::unique_ptr<tKeyType[]> tArrayPtr;

private:
    std::size_t m_arraySize;
    std::size_t m_capacity;

    tArrayPtr m_array;

    std::size_t m_ones;

private:
    static const uint cENTRIES = (CHAR_BIT * sizeof(tKeyType));
};

class Concurrent_Bit_Set {

    friend class Bit_Set;

public:
    typedef uint tKeyType;

public:

    Concurrent_Bit_Set(std::size_t size)
            : m_array(nullptr),
              m_ones(0) {

        // Initialize cells
        m_arraySize = (size + cENTRIES - 1) / cENTRIES; //get the ceiling
        m_capacity = m_arraySize * cENTRIES;
        try {
            m_array.reset(new std::atomic<tKeyType>[m_arraySize]()); //value init
        } catch (std::bad_alloc &e) {
            std::cerr << e.what() << std::endl;
            raise(SIGINT);
        }
    }

    Concurrent_Bit_Set(Concurrent_Bit_Set &&other)
            : m_arraySize(other.m_arraySize),
              m_capacity(other.m_capacity),
              m_array(std::move(other.m_array)),
              m_ones(other.m_ones.load()) { }

    //Conversion
    Concurrent_Bit_Set(Bit_Set &&other);

    Concurrent_Bit_Set &operator=(Concurrent_Bit_Set &&other) {
        m_arraySize = other.m_arraySize;
        m_capacity = other.m_capacity;
        m_array = std::move(other.m_array);
        m_ones.store(other.m_ones.load());

        return *this;
    }

    bool testAndSet(const tKeyType &idx) {
        RAISE(idx / cENTRIES < m_arraySize);
        RAISE(idx < m_capacity);

        std::size_t i = idx / cENTRIES;
        tKeyType mask = 1 << (idx % cENTRIES);

        tKeyType des, val = m_array[i].load();

        while (!(val & mask)) {
            des = val | mask;

            bool cas = m_array[i].compare_exchange_strong(val, des);

            if (cas) {
                ++m_ones; // we just set the bit, inc m_ones
                return false; // the bit was _not_ set before
            }
        }

        // value was already set
        return true;
    }

    bool isSet(const tKeyType &idx) const {
        RAISE(idx / cENTRIES < m_arraySize);
        RAISE(idx < m_capacity);

        tKeyType val = m_array[idx / cENTRIES].load();
        return val && (val & (1 << (idx % cENTRIES)));
    }

    bool zero() const { return m_ones.load() == 0; }

    bool one() const { return m_ones.load() == m_capacity; }

    bool zero(const std::size_t idx) const {
        return !isSet(idx);
    }

    bool one(const std::size_t idx) const {
        return isSet(idx);
    }

    auto capacity() const { return m_capacity; }

    auto zeros() const { return m_capacity - m_ones.load(); }

    auto ones() const { return m_ones.load(); }

    template<class BitSet>
    void unsafe_copy(const BitSet &other) {

        m_arraySize = other.m_arraySize;
        m_capacity = other.m_capacity;
        m_ones.store(other.m_ones);

        try {
            m_array.reset(new std::atomic<tKeyType>[m_arraySize]); //random init
        } catch (std::bad_alloc &e) {
            std::cerr << e.what() << std::endl;
            raise(SIGINT);
        }

        for (std::size_t i = 0; i < m_arraySize; ++i) {
            m_array[i].store(other.m_array[i]);
        }

    }

    template<class BitSet>
    void unsafe_merge(BitSet &&other) {

        std::size_t bSize;
        tArrayPtr bArray;
        if (m_capacity >= other.m_capacity) {
            bSize = other.m_arraySize;
            bArray = std::move(other.m_array);
        } else {
            bSize = m_arraySize;
            bArray = std::move(m_array);

            m_arraySize = other.m_arraySize;
            m_capacity = other.m_capacity;
            m_array = std::move(other.m_array);
        }

        m_ones.store(0);
        for (std::size_t i = 0; i < bSize; ++i) {
            m_array[i].fetch_or(bArray[i].load());
            m_ones.fetch_add(popcount(m_array[i].load()));
        }

    }

    template<class BitSet>
    void unsafe_filter(const BitSet &filter) {

        m_ones.store(0);
        for (std::size_t i = 0; i < std::min(m_arraySize, filter.m_arraySize); ++i) {
            m_array[i].fetch_and(~filter.m_array[i].load());
            m_ones.fetch_add(popcount(m_array[i].load()));
        }
    }

    // compat definitions
    bool insert(const tKeyType &key) { return !testAndSet(key); }

    bool contains(const tKeyType &key) const { return isSet(key); }

    std::size_t count(const tKeyType &key) const { return contains(key); }

    bool empty() const { return zero(); };

    bool empty(const std::size_t idx) const { return zero(idx); }

    tKeyType at(const std::size_t idx) const { return idx * isSet(idx); }

    auto size() const { return ones(); }


public:
    typedef _detail::iterator<Concurrent_Bit_Set, tKeyType> iterator;
    typedef _detail::range_type<Concurrent_Bit_Set, iterator> range_type;

    iterator begin() const {
        return iterator(*this, 0, true);
    }

    iterator end() const {
        return iterator(*this, m_capacity, false);
    }

    range_type range() const {
        return range_type(*this);
    }

private:
    typedef std::unique_ptr<std::atomic<tKeyType>[]> tArrayPtr;

private:
    std::size_t m_arraySize;
    std::size_t m_capacity;

    tArrayPtr m_array;

    std::atomic<std::size_t> m_ones;

private:
    static const uint cENTRIES = (CHAR_BIT * sizeof(tKeyType));
};