#pragma once

#include <atomic>
#include <mutex>
#include <condition_variable>
#include <functional>
#include <cstring>
#include <csignal>

#include <iostream>
#include <iomanip>

#include "utils/Misc.h"
#include "utils/ASSERT.h"
#include "utils/Random.h"
#include "utils/VTuneAdapter.h"

#include <tbb/tbb_stddef.h>
#include <tbb/parallel_for.h>
#include <tbb/blocked_range.h>

#include <boost/serialization/access.hpp>
#include <boost/serialization/array.hpp>

template<typename T>
struct Hasher {

    Hasher() {
        std::uniform_int_distribution<T> dist(std::numeric_limits<T>::min(), std::numeric_limits<T>::max());
        while (a = dist(startGen), a % 2 == 0) { } // find odd value
    }

    T operator()(const T x) const {
        return (a * x);
    }

    T a;
};

struct InsertReturn {
public:
    enum class State : unsigned char {
        False, True, Full
    };

    InsertReturn() : m_state(State::False) { }

    InsertReturn(bool b) : m_state(b ? State::True : State::False) { }

    InsertReturn(const InsertReturn::State &s) : m_state(s) { }

    InsertReturn(const InsertReturn &o) : m_state(o.m_state) { }

    operator bool() const {
        return m_state == State::True;
    }

    bool operator==(const InsertReturn::State &other) const { return m_state == other; }

private:
    State m_state;
};

namespace _detail {

    template<class Container, typename Value>
    struct iterator : public std::iterator<std::bidirectional_iterator_tag, Value> {

    public:
        iterator(const Container &container, std::size_t idx, bool findNext)
                : m_container(container),
                  m_idx(idx) {

            if (findNext)
               _advance(true);

        }

        iterator(const iterator &o)
                : m_container(o.m_container),
                  m_idx(o.m_idx) { }

        iterator &operator++() {
            _advance(false);
            return *this;
        }

        iterator operator++(int) {
            auto tmp = *this;
            ++(*this);
            return tmp;
        }

        iterator &operator--() {
            _dec(false);
            return *this;
        }

        iterator operator--(int) {
            auto tmp = *this;
            --(*this);
            return tmp;
        }

        std::size_t operator-(const iterator &other) const {
            return m_idx - other.m_idx;
        }

        std::size_t operator+(const iterator &other) const {
            return m_idx + other.m_idx;
        }

        iterator operator=(const iterator &other) {
            m_idx = other.m_idx;
            return *this;
        }

        void setIdx(const std::size_t &idx, bool findNext) {
            m_idx = idx;
            if (findNext)
                _advance(true);
        }

        bool operator==(const iterator &j) const {
            return m_idx == j.m_idx;
        }

        bool operator!=(const iterator &j) const { return !(*this == j); }

        bool operator<(const iterator &other) const {
            return m_idx < other.m_idx;
        }

        bool operator<=(const iterator &other) const {
            return m_idx <= other.m_idx;
        }

        bool operator>(const iterator &other) const {
            return m_idx > other.m_idx;
        }

        bool operator>=(const iterator &other) const {
            return m_idx >= other.m_idx;
        }

        auto operator*() const { return m_container.at(m_idx); }

        //const auto *operator->() const { return &m_container.at(m_idx); }

    private:
        void _advance(const bool testFirst){
            if(testFirst){
                while ((m_idx < m_container.capacity() && m_container.empty(m_idx)))
                    ++m_idx;
            } else {
                while (++m_idx, (m_idx < m_container.capacity() && m_container.empty(m_idx))) { }
            }
        }

        void _dec(const bool testFirst){
            if(testFirst){
                while ((m_idx >= 0 && m_container.empty(m_idx)))
                    --m_idx;
            } else {
                while (--m_idx, (m_idx >= 0 && m_container.empty(m_idx))) { }
            }
        }

    protected:
        const Container &m_container;
        std::size_t m_idx;
    };

    template<class Container, class IT>
    class range_type {
        const Container &m_container;
        IT m_begin;
        IT m_end;
        mutable IT m_midpoint;
        std::size_t m_grainsize;

    public:

        //! True if range is empty.
        bool empty() const { return m_end <= m_begin; }

        //! True if range can be partitioned into two subranges.
        bool is_divisible() const {
            return m_end - m_begin >= grainsize();
        }

        //! Split range.
        range_type(range_type &r, tbb::split) :
                m_container(r.m_container),
                m_begin(r.m_begin),
                m_end(r.m_midpoint),
                m_midpoint(m_begin),
                m_grainsize(r.m_grainsize){

            r.m_begin = r.m_midpoint;

            set_midpoint();
            r.set_midpoint();
        }

        //! Init range with container and grainsize specified
        range_type(const Container &cont, const std::size_t & grainsize = 1e3) :
                m_container(cont),
                m_begin(cont.begin()),
                m_end(cont.end()),
                m_midpoint(cont.begin()),
                m_grainsize(grainsize){
            set_midpoint();
        }

        const IT &begin() const { return m_begin; }

        const IT &end() const { return m_end; }

        //! The grain size for this range.
        std::size_t grainsize() const { return m_grainsize; }

        //! Set my_midpoint_node to point approximately half way between my_begin_node and my_end_node.
        void set_midpoint() const {
            if (m_begin == m_end) // not divisible
                m_midpoint = m_end;
            else {
                m_midpoint.setIdx((m_end + m_begin) / 2, true);
            }
        }
    };
}

//forward declare concurrent version
class Concurrent_LP_Set;

class LP_Set {

    friend class Concurrent_LP_Set;

public:
    typedef uint tKeyType;

public:

    LP_Set(std::size_t size = 1, Hasher<tKeyType> hasher = Hasher<tKeyType>())
            : m_items(0),
              m_array(nullptr),
              m_hasher(hasher) {
        // Initialize cells
        m_arraySize = nextPow2(size);
        try {
            m_array.reset(new tKeyType[m_arraySize]()); //value init
        } catch (std::bad_alloc &e) {
            std::cerr << e.what() << std::endl;
            raise(SIGINT);
        }
    }

    LP_Set(LP_Set &&other)
            : m_arraySize(other.m_arraySize),
              m_items(other.m_items),
              m_array(std::move(other.m_array)),
              m_hasher(std::move(other.m_hasher)) { }

    //Conversion
    LP_Set(Concurrent_LP_Set &&other);

    LP_Set &operator=(LP_Set &&other) {
        m_arraySize = other.m_arraySize;
        m_items = other.m_items;
        m_array = std::move(other.m_array);
        m_hasher = std::move(other.m_hasher);

        return *this;
    }

    //Copy function
    LP_Set copy() const {
        LP_Set c(m_arraySize, m_hasher);
        c.m_items = m_items;

        for(uint i = 0; i < m_arraySize; ++i){
            c.m_array[i] = m_array[i];
        }

        return c;
    }

    bool insert(const tKeyType &key) {
        ASSERT(key != 0);

        if (!m_rehashing && (double) m_items / (double) m_arraySize > 0.5)
            rehash(m_arraySize << 1);

        for (tKeyType idx = m_hasher(key); ; idx++) {
            idx &= m_arraySize - 1;
            ASSERT(idx < m_arraySize);

            // Load the key that was there.
            tKeyType probedKey = m_array[idx];

            if (probedKey == key)
                return false; // the key is already in the set, return false;
            else {
                // The entry was either free, or contains another key.
                if (probedKey != 0)
                    continue; // Usually, it contains another key. Keep probing.
                // The entry was free. take it
                m_array[idx] = key;
                ++m_items;
                return true;
            }
        }
    }

    bool contains(const tKeyType &key) const {
        ASSERT(key != 0);

        for (tKeyType idx = m_hasher(key); ; idx++) {
            idx &= m_arraySize - 1;
            tKeyType probedKey = m_array[idx];
            if (probedKey == key)
                return true;;
            if (probedKey == 0)
                return false;
        }

    }

    std::size_t count(const tKeyType &key) const {
        return contains(key);
    }

    bool empty() const { return m_items == 0; };

    bool empty(const std::size_t idx) const {
        return m_array[idx] == 0;
    };

    tKeyType at(const std::size_t idx) const {

        return m_array[idx];
    }

    auto capacity() const { return m_arraySize; }

    auto size() const { return m_items; }

    void rehash(std::size_t newSize) {

        m_rehashing = true;

        std::size_t oldSize = m_arraySize;
        m_arraySize = nextPow2(newSize);

        auto oldArray = std::move(m_array);
        m_items = 0;

        try {
            m_array.reset(new tKeyType[m_arraySize]()); //value init
        } catch (std::bad_alloc &e) {
            std::cerr << e.what() << std::endl;
            raise(SIGINT);
        }

        for (std::size_t i = 0; i < oldSize; ++i) {
            if (oldArray[i] != 0)
                insert(oldArray[i]);
        }

        m_rehashing = false;

    }

    void merge(LP_Set &&other) {

        rehash((capacity() + other.capacity()));

        for (std::size_t i = 0; i < other.capacity(); ++i) {
            if (other.m_array[i] != 0)
                insert(other.m_array[i]);
        }
    }

    template<class Set>
    void rehash(std::size_t newSize, const Set &filter) {

        m_rehashing = true;

        std::size_t oldSize = m_arraySize;
        m_arraySize = nextPow2(newSize);

        auto oldArray = std::move(m_array);
        m_items = 0;

        try {
            m_array.reset(new tKeyType[m_arraySize]()); //value init
        } catch (std::bad_alloc &e) {
            std::cerr << e.what() << std::endl;
            raise(SIGINT);
        }

        for (std::size_t i = 0; i < oldSize; ++i) {
            if (oldArray[i] != 0 && !filter.count(oldArray[i]))
                insert(oldArray[i]);
        }

        m_rehashing = false;
    }

    template<class Set>
    void merge(LP_Set &&other, const Set &filter) {
        rehash((capacity() + other.capacity()), filter);

        for (std::size_t i = 0; i < other.capacity(); ++i) {
            if (other.m_array[i] != 0 && !filter.count(other.m_array[i]))
                insert(other.m_array[i]);
        }
    }


public:
    typedef _detail::iterator<LP_Set, tKeyType> iterator;
    typedef _detail::range_type<LP_Set, iterator> range_type;

    iterator begin() const {
        return iterator(*this, 0, true);
    }

    iterator end() const {
        return iterator(*this, m_arraySize, false);
    }

    range_type range() const {
        return range_type(*this);
    }


private:
    std::size_t m_arraySize;
    std::size_t m_items;
    std::unique_ptr<tKeyType[]> m_array;
    Hasher<tKeyType> m_hasher;

    bool m_rehashing = false;

private:
    friend class boost::serialization::access;

    template<class Archive>
    void save(Archive & ar, __attribute__((unused)) const unsigned int version) const
    {
        // invoke serialization of the base class
        ar << m_arraySize;
        ar << m_items;
        ar << boost::serialization::make_array<tKeyType>(m_array.get(), m_arraySize);
        ar << m_hasher.a;
    }

    template<class Archive>
    void load(Archive & ar, __attribute__((unused)) const unsigned int version)
    {
        // invoke serialization of the base class
        ar >> m_arraySize;
        ar >> m_items;

        try {
            m_array.reset(new tKeyType[m_arraySize]); //random init
        } catch (std::bad_alloc &e) {
            std::cerr << e.what() << std::endl;
            raise(SIGINT);
        }
        ar >> boost::serialization::make_array<tKeyType>(m_array.get(), m_arraySize);
        ar >> m_hasher.a;
    }

    BOOST_SERIALIZATION_SPLIT_MEMBER()

public:
    template <class Set>
    bool operator==(const Set & set) const {
        if(size() != set.size())
            return false;

        for(const auto & i : set){
            if(!contains(i))
                return false;
        }

        return true;
    }

};

template<class HT>
class GrowingHashTableHandle;

template<class HT>
class ConstGrowingHashTableHandle;


template<class HT>
class GrowingHashTable;

class Concurrent_LP_Set {

    friend class LP_Set;

    friend class GrowingHashTable<Concurrent_LP_Set>;

    friend class GrowingHashTableHandle<Concurrent_LP_Set>;

    friend class ConstGrowingHashTableHandle<Concurrent_LP_Set>;

public:
    typedef uint tKeyType;

public:

    Concurrent_LP_Set(const std::size_t size, const uint version = 0,
                      Hasher<tKeyType> hasher = Hasher<tKeyType>())
            : m_items(0),
              m_array(nullptr),
              m_version(version),
              m_hasher(hasher),
              m_currentCopyBlock(0) {
        // Initialize cells
        m_arraySize = nextPow2(size);
        try {
            m_array.reset(new std::atomic<tKeyType>[m_arraySize]()); //value init
        } catch (std::bad_alloc &e) {
            std::cerr << e.what() << std::endl;
            raise(SIGINT);
        }
    }

    Concurrent_LP_Set(Concurrent_LP_Set &&other)
            : m_arraySize(other.m_arraySize),
              m_items(other.m_items.load()),
              m_array(std::move(other.m_array)),
              m_version(other.m_version),
              m_hasher(std::move(other.m_hasher)),
              m_currentCopyBlock(0) { }

    //conversion
    Concurrent_LP_Set(LP_Set &&other);

    Concurrent_LP_Set &operator=(Concurrent_LP_Set &&other) {
        m_arraySize = other.m_arraySize;
        m_items.store(other.m_items.load());
        m_array = std::move(other.m_array);
        m_version = other.m_version;
        m_hasher = std::move(other.m_hasher);
        m_currentCopyBlock.store(other.m_currentCopyBlock.load());

        return *this;
    }


    InsertReturn insert(const tKeyType &key) {
        ASSERT(key != 0);

        if (m_items.load() > m_arraySize >> 1)
            return InsertReturn::State::Full;

        uint cols = 0;
        for (tKeyType idx = m_hasher(key); ; idx++) {

            if (cols > c_growThreshold) {
                return InsertReturn::State::Full;
            }

            idx &= m_arraySize - 1;
            ASSERT(idx < m_arraySize);

            // Load the key that was there.
            tKeyType probedKey = m_array[idx].load();

            if (probedKey == key) {
                return false; // the key is already in the set, return false;
            }
            else {
                // The entry was either free, or contains another key.
                if (probedKey != 0) {
                    ++cols;
                    continue; // Usually, it contains another key. Keep probing.
                }
                // The entry was free. Now let's try to take it using a CAS.
                tKeyType prevKey = 0;
                bool cas = m_array[idx].compare_exchange_strong(prevKey, key);
                if (cas) {
                    ++m_items;
                    return true; // we just added the key to the set
                }
                else if (prevKey == key) {
                    return false; // the key was already added by another thread
                }
                else {
                    ++cols;
                    continue; // another thread inserted a different key in this position
                }
            }
        }
    }

    tKeyType get(const tKeyType &key) const {
        ASSERT(key != 0);

        for (tKeyType idx = m_hasher(key); ; idx++) {
            idx &= m_arraySize - 1;
            tKeyType probedKey = m_array[idx].load();
            if (probedKey == key)
                return key;;
            if (probedKey == 0)
                return 0;
        }

    }

    bool contains(const tKeyType &key) const {
        return get(key);
    }

    std::size_t count(const tKeyType &key) const {
        return contains(key);
    }

    bool empty() const { return m_items.load() == 0; };

    bool empty(const std::size_t idx) const { return m_array[idx].load(std::memory_order_relaxed) == 0; };

    tKeyType at(const std::size_t idx) const { return m_array[idx].load(std::memory_order_relaxed); }

    auto capacity() const { return m_arraySize; }

    auto size() const { return m_items.load(); }

    template<class Source>
    void unsafe_copy(const Source &other) {
        ASSERT(m_arraySize == other.capacity());
        for (std::size_t i = 0; i < m_arraySize; ++i)
            m_array[i] = other.at(i);
        m_items.store(other.size());
    }

    template<class Source, class Filter>
    void unsafe_merge(Source &&other, const Filter &filter) {

        VTUNE_TASK(SetMergeAllocate);
        std::size_t oldSize = m_arraySize;

        std::size_t combinedItems = size() + other.size();
        std::size_t newSize = nextPow2(combinedItems);

        while (combinedItems > newSize >> 1)
            newSize <<= 1;

        m_arraySize = newSize;

        auto oldArray = std::move(m_array);
        m_items.store(0, std::memory_order_relaxed);

        try {
            m_array.reset(new std::atomic<tKeyType>[m_arraySize]); //random init
        } catch (std::bad_alloc &e) {
            std::cerr << e.what() << std::endl;
            raise(SIGINT);
        }
        VTUNE_END_TASK(SetMergeAllocate);

        VTUNE_TASK(SetMergeFill);
        tbb::parallel_for(tbb::blocked_range<std::size_t>(0, m_arraySize), [this](const auto &r) {
            for (auto i = r.begin(); i != r.end(); ++i) {
                m_array[i].store(0, std::memory_order_relaxed);
            }
        });
        VTUNE_END_TASK(SetMergeFill);

        VTUNE_TASK(SetMergeCopy);
        tbb::parallel_for(tbb::blocked_range<std::size_t>(0, std::max(oldSize, other.capacity())),
                          [&oldArray, oldSize, &other, &filter, this](const auto &r) {

                              tKeyType val = 0;
                              for (auto i = r.begin(); i != r.end(); ++i) {
                                  if (i < oldSize) {
                                      val = oldArray[i].load(std::memory_order_relaxed);

                                      if (val != 0 && !filter.count(val))
                                          this->_insert(val);
                                  }

                                  if (i < other.capacity()) {
                                      val = other.at(i);

                                      if (val != 0 && !filter.count(val))
                                          this->_insert(val);
                                  }
                              }
                          });
    }

private:

    void _insert(const tKeyType &key) {

        ASSERT(key != 0);

        for (tKeyType idx = m_hasher(key); ; idx++) {
            idx &= m_arraySize - 1;
            ASSERT(idx < m_arraySize);

            // Load the key that was there.
            tKeyType probedKey = m_array[idx].load();

            if (probedKey == key) {
                return; // the key is already in the set, return false;
            }
            else {
                // The entry was either free, or contains another key.
                if (probedKey != 0)
                    continue; // Usually, it contains another key. Keep probing.
                // The entry was free. Now let's try to take it using a CAS.
                tKeyType prevKey = 0;
                bool cas = m_array[idx].compare_exchange_strong(prevKey, key);
                if (cas) {
                    ++m_items;
                    return; // we just added the key to the set
                }
                else if (prevKey == key) {
                    return; // the key was already added by another thread
                }
                else
                    continue; // another thread inserted a different key in this position
            }
        }
    }

    void _migrate(const std::size_t idx, Concurrent_LP_Set &target) {
        auto val = at(idx);
        if (val != 0)
            target._insert(val);
    }

public:
    typedef _detail::iterator<Concurrent_LP_Set, tKeyType> iterator;
    typedef _detail::range_type<Concurrent_LP_Set, iterator> range_type;

    iterator begin() const {
        return iterator(*this, 0, true);
    }

    iterator end() const {
        return iterator(*this, m_arraySize, false);
    }

    range_type range() const {
        return range_type(*this);
    }


private:
    std::size_t m_arraySize;
    std::atomic<std::size_t> m_items;
    std::unique_ptr<std::atomic<tKeyType>[]> m_array;
    uint m_version;
    Hasher<tKeyType> m_hasher;

    std::atomic<std::size_t> m_currentCopyBlock;

    const static uint c_growThreshold = 10;

};