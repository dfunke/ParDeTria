#pragma once

#include <atomic>
#include <functional>
#include <cstring>
#include <csignal>

#include <iostream>
#include <iomanip>

#include "Misc.h"
#include "ASSERT.h"
#include "Random.h"

#include <tbb/tbb_stddef.h>
#include <tbb/parallel_for.h>
#include <tbb/blocked_range.h>

template<typename T>
struct Hasher {

    Hasher() {
        std::uniform_int_distribution<T> dist(std::numeric_limits<T>::min(), std::numeric_limits<T>::max());
        while (a = dist(startGen), a % 2 == 0) { } // find odd value

        l = sizeof(T) * CHAR_BIT;
    }

    T operator()(const T x) const {
        return (a * x) >> (sizeof(T) * CHAR_BIT - l);
    }

    T a;
    T l;
};

namespace _detail {

    template<class Container, typename Value>
    struct iterator : public std::iterator<std::bidirectional_iterator_tag, Value> {

    public:
        iterator(const Container &container, std::size_t idx, bool findNext)
                : m_container(container),
                  m_idx(idx) {

            if (findNext)
                while ((m_idx < m_container.capacity() && m_container.empty(m_idx)))
                    ++m_idx;

        }

        iterator(const iterator &o)
                : m_container(o.m_container),
                  m_idx(o.m_idx) { }

        iterator &operator++() {
            while (++m_idx, (m_idx < m_container.capacity() && m_container.empty(m_idx))) { }
            return *this;
        }

        iterator operator++(int) {
            auto tmp = *this;
            ++(*this);
            return tmp;
        }

        iterator &operator--() {
            while (--m_idx, (m_idx >= 0 && m_container.empty(m_idx))) { }
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
                while ((m_idx < m_container.capacity() && m_container.empty(m_idx)))
                    ++m_idx;
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

        const auto operator*() const { return m_container.at(m_idx); }

        //const auto *operator->() const { return &m_container.at(m_idx); }

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
                m_midpoint(m_begin) {

            r.m_begin = r.m_midpoint;

            set_midpoint();
            r.set_midpoint();
        }

        //! Init range with container and grainsize specified
        range_type(const Container &cont) :
                m_container(cont),
                m_begin(cont.begin()),
                m_end(cont.end()),
                m_midpoint(cont.begin()) {
            set_midpoint();
        }

        const IT &begin() const { return m_begin; }

        const IT &end() const { return m_end; }

        //! The grain size for this range.
        std::size_t grainsize() const { return 1; }

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

class LP_Set {

public:
    typedef uint tKeyType;

public:

    LP_Set(std::size_t size, Hasher<tKeyType> hasher = Hasher<tKeyType>())
            : m_items(0),
              m_array(nullptr),
              m_hasher(hasher) {
        // Initialize cells
        m_arraySize = nextPow2(size);
        m_array = std::make_unique<std::vector<tKeyType>>();
        m_array->resize(m_arraySize);

        m_hasher.l = log2(m_arraySize);
    }

    LP_Set(LP_Set &&other)
            : m_arraySize(other.m_arraySize),
              m_items(other.m_items),
              m_array(std::move(other.m_array)),
              m_hasher(std::move(other.m_hasher)) { }

    LP_Set &operator=(LP_Set &&other) {
        m_arraySize = other.m_arraySize;
        m_items = other.m_items;
        m_array = std::move(other.m_array);
        m_hasher = std::move(other.m_hasher);

        return *this;
    }

    ~LP_Set() { }

    bool insert(const tKeyType &key) {
        ASSERT(key != 0);

        if (!m_rehashing && m_items / m_arraySize > 0.5)
            rehash(m_arraySize << 1);

        for (tKeyType idx = m_hasher(key); ; idx++) {
            idx &= m_arraySize - 1;
            ASSERT(idx < m_arraySize);

            // Load the key that was there.
            tKeyType probedKey = m_array->at(idx);

            if (probedKey == key)
                return false; // the key is already in the set, return false;
            else {
                // The entry was either free, or contains another key.
                if (probedKey != 0)
                    continue; // Usually, it contains another key. Keep probing.
                // The entry was free. take it
                m_array->at(idx) = key;
                ++m_items;
                return true;
            }
        }
    }

    bool contains(const tKeyType &key) const {
        ASSERT(key != 0);

        for (tKeyType idx = m_hasher(key); ; idx++) {
            idx &= m_arraySize - 1;
            tKeyType probedKey = m_array->at(idx);
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
        return m_array->at(idx) == 0;
    };

    tKeyType at(const std::size_t idx) const {

        return m_array->at(idx);
    }

    auto capacity() const { return m_arraySize; }

    auto size() const { return m_items; }

    void rehash(std::size_t newSize) {

        m_rehashing = true;

        std::size_t oldSize = m_arraySize;
        m_arraySize = nextPow2(newSize);
        m_hasher.l = log2(m_arraySize);

        auto oldArray = std::move(m_array);
        m_items = 0;

        m_array = std::make_unique<std::vector<tKeyType>>();
        m_array->resize(m_arraySize);

        for (std::size_t i = 0; i < oldSize; ++i) {
            if (oldArray->at(i) != 0)
                insert(oldArray->at(i));
        }

        m_rehashing = false;

    }

    void merge(LP_Set &&other) {

        rehash((capacity() + other.capacity()) << 1);

        for (std::size_t i = 0; i < other.capacity(); ++i) {
            if (other.m_array->at(i) != 0)
                insert(other.m_array->at(i));
        }
    }

    template<class Set>
    void rehash(std::size_t newSize, const Set &filter) {

        m_rehashing = true;

        std::size_t oldSize = m_arraySize;
        m_arraySize = nextPow2(newSize);
        m_hasher.l = log2(m_arraySize);

        auto oldArray = std::move(m_array);
        m_items = 0;

        m_array = std::make_unique<std::vector<tKeyType>>();
        m_array->resize(m_arraySize);

        for (std::size_t i = 0; i < oldSize; ++i) {
            if (oldArray->at(i) != 0 && !filter.count(oldArray->at(i)))
                insert(oldArray->at(i));
        }

        m_rehashing = false;
    }

    template<class Set>
    void merge(LP_Set &&other, const Set &filter) {
        rehash((capacity() + other.capacity()) << 1, filter);

        for (std::size_t i = 0; i < other.capacity(); ++i) {
            if (other.m_array->at(i) != 0 && !filter.count(other.m_array->at(i)))
                insert(other.m_array->at(i));
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
    std::unique_ptr<std::vector<tKeyType>> m_array;
    Hasher<tKeyType> m_hasher;

    bool m_rehashing = false;

};

class Concurrent_LP_Set {

public:
    typedef uint tKeyType;

public:

    Concurrent_LP_Set(std::size_t size, Hasher<tKeyType> hasher = Hasher<tKeyType>())
            : m_items(0),
              m_array(nullptr),
              m_hasher(hasher) {
        // Initialize cells
        m_arraySize = nextPow2(size);
        m_array = std::unique_ptr<std::atomic<tKeyType>[]>(new std::atomic<tKeyType>[m_arraySize]());

        m_hasher.l = log2(m_arraySize);
    }

    Concurrent_LP_Set(Concurrent_LP_Set &&other)
            : m_arraySize(other.m_arraySize),
              m_items(other.m_items.load(std::memory_order_relaxed)),
              m_array(std::move(other.m_array)),
              m_hasher(std::move(other.m_hasher)) { }

    Concurrent_LP_Set &operator=(Concurrent_LP_Set &&other) {
        m_arraySize = other.m_arraySize;
        m_items.store(other.m_items.load(std::memory_order_relaxed), std::memory_order_relaxed);
        m_array = std::move(other.m_array);
        m_hasher = std::move(other.m_hasher);

        return *this;
    }

    ~Concurrent_LP_Set() { }

    bool insert(const tKeyType &key) {
        ASSERT(key != 0);

        if (m_items.load(std::memory_order_relaxed) / m_arraySize > 0.75)
            throw std::length_error("Overfull Concurrent Set");

        for (tKeyType idx = m_hasher(key); ; idx++) {
            idx &= m_arraySize - 1;
            ASSERT(idx < m_arraySize);

            // Load the key that was there.
            tKeyType probedKey = m_array[idx].load(std::memory_order_relaxed);

            if (probedKey == key)
                return false; // the key is already in the set, return false;
            else {
                // The entry was either free, or contains another key.
                if (probedKey != 0)
                    continue; // Usually, it contains another key. Keep probing.
                // The entry was free. Now let's try to take it using a CAS.
                tKeyType prevKey = 0;
                bool cas = m_array[idx].compare_exchange_strong(prevKey, key, std::memory_order_relaxed);
                if (cas) {
                    ++m_items;
                    return true; // we just added the key to the set
                }
                else if (prevKey == key)
                    return false; // the key was already added by another thread
                else
                    continue; // another thread inserted a different key in this position
            }
        }
    }

    bool contains(const tKeyType &key) const {
        ASSERT(key != 0);
        for (tKeyType idx = m_hasher(key); ; idx++) {
            idx &= m_arraySize - 1;
            tKeyType probedKey = m_array[idx].load(std::memory_order_relaxed);
            if (probedKey == key)
                return true;;
            if (probedKey == 0)
                return false;
        }

    }

    std::size_t count(const tKeyType &key) const {
        return contains(key);
    }

    bool empty() const { return m_items.load(std::memory_order_relaxed) == 0; };

    bool empty(const std::size_t idx) const { return m_array[idx].load(std::memory_order_relaxed) == 0; };

    tKeyType at(const std::size_t idx) const { return m_array[idx].load(std::memory_order_relaxed); }

    auto capacity() const { return m_arraySize; }

    auto size() const { return m_items.load(); }

    void unsafe_rehash(std::size_t newSize) {
        std::size_t oldSize = m_arraySize;
        m_arraySize = nextPow2(newSize);
        m_hasher.l = log2(m_arraySize);

        auto oldArray = std::move(m_array);
        m_items.store(0, std::memory_order_relaxed);

        m_array = std::unique_ptr<std::atomic<tKeyType>[]>(new std::atomic<tKeyType>[m_arraySize]());

        tbb::parallel_for(std::size_t(0), oldSize, [&oldArray, this](const uint i) {
            if (oldArray[i].load(std::memory_order_relaxed) != 0)
                insert(oldArray[i].load(std::memory_order_relaxed));
        });
    }

    void unsafe_merge(Concurrent_LP_Set &&other) {

        unsafe_rehash((capacity() + other.capacity()) << 1);

        tbb::parallel_for(std::size_t(0), other.capacity(), [&other, this](const uint i) {
            if (other.m_array[i].load(std::memory_order_relaxed) != 0)
                insert(other.m_array[i].load(std::memory_order_relaxed));
        });
    }


    template<class Set>
    void unsafe_merge(Concurrent_LP_Set &&other, const Set &filter) {

        std::size_t oldSize = m_arraySize;
        m_arraySize = nextPow2((capacity() + other.capacity()) << 1);
        m_hasher.l = log2(m_arraySize);

        auto oldArray = std::move(m_array);
        m_items.store(0, std::memory_order_relaxed);

        m_array = std::unique_ptr<std::atomic<tKeyType>[]>(new std::atomic<tKeyType>[m_arraySize]());

        tbb::parallel_for(tbb::blocked_range<std::size_t>(0, std::max(oldSize, other.capacity())), [&oldArray, oldSize, &other, &filter, this](const auto & r) {

            tKeyType  val = 0;
            for(auto i = r.begin(); i != r.end(); ++i) {
                if (i < oldSize){
                    val = oldArray[i].load(std::memory_order_relaxed);

                    if(val != 0 && !filter.count(val))
                        this->insert(val);
                }

                if (i < other.capacity()){
                    val = other.m_array[i].load(std::memory_order_relaxed);

                    if(val != 0 && !filter.count(val))
                        this->insert(val);
                }
            }
        });
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
    Hasher<tKeyType> m_hasher;

};
