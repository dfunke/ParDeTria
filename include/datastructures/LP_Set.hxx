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

        auto operator*() const { return m_container.at(m_idx); }

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

        rehash((capacity() + other.capacity()));

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
        rehash((capacity() + other.capacity()), filter);

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
              m_hasher(hasher),
              m_fGrowing(false) {
        // Initialize cells
        m_arraySize = nextPow2(size);
        try {
            m_array.reset(new std::atomic<tKeyType>[m_arraySize]()); //value init
        } catch (std::bad_alloc &e) {
            std::cerr << e.what() << std::endl;
            raise(SIGINT);
        }

        m_hasher.l = log2(m_arraySize);
    }

    Concurrent_LP_Set(Concurrent_LP_Set &&other)
            : m_arraySize(other.m_arraySize),
              m_items(other.m_items.load()),
              m_array(std::move(other.m_array)),
              m_hasher(std::move(other.m_hasher)),
              m_fGrowing(false) { }

    Concurrent_LP_Set &operator=(Concurrent_LP_Set &&other) {
        m_arraySize = other.m_arraySize;
        m_items.store(other.m_items.load());
        m_array = std::move(other.m_array);
        m_hasher = std::move(other.m_hasher);

        return *this;
    }


    bool insert(const tKeyType &key) {
        ASSERT(key != 0);

        if (m_fGrowing.load())
            helpGrowing();

        bool result = false;
        std::size_t steps = 0; // count number of steps
        auto currArr = m_array.get(); // store current array pointer for later comparision

        for (tKeyType idx = m_hasher(key); ; idx++) {

            if (steps > m_arraySize / 4) {
                //we stepped through more than a quarter the array > grow
                grow();

                //reset insertion
                steps = 0;
                idx = m_hasher(key);
            }

            idx &= m_arraySize - 1;
            ASSERT(idx < m_arraySize);

            // Load the key that was there.
            tKeyType probedKey = m_array[idx].load();

            if (probedKey == key) {
                result = false; // the key is already in the set, return false;
                break;
            }
            else {
                // The entry was either free, or contains another key.
                if (probedKey != 0) {
                    ++steps;
                    continue; // Usually, it contains another key. Keep probing.
                }
                // The entry was free. Now let's try to take it using a CAS.
                tKeyType prevKey = 0;
                bool cas = m_array[idx].compare_exchange_strong(prevKey, key);
                if (cas) {
                    ++m_items;
                    result = true; // we just added the key to the set
                    break;
                }
                else if (prevKey == key) {
                    result = false; // the key was already added by another thread
                    break;
                }
                else {
                    ++steps;
                    continue; // another thread inserted a different key in this position
                }
            }
        }

        if (result && (m_fGrowing.load() || currArr != m_array.get())) {
            //we inserted a key but the array has changed or we are currently growing -> re-insert element
            //TODO ABA problem

            helpGrowing();
            migrate(key, m_array, m_arraySize, m_hasher);
        }

        return result;
    }

    bool contains(const tKeyType &key) const {
        ASSERT(key != 0);

        if (m_fGrowing.load())
            helpGrowing();

        for (tKeyType idx = m_hasher(key); ; idx++) {
            idx &= m_arraySize - 1;
            tKeyType probedKey = m_array[idx].load();
            if (probedKey == key)
                return true;;
            if (probedKey == 0)
                return false;
        }

    }

    std::size_t count(const tKeyType &key) const {
        return contains(key);
    }

    bool empty() const { return m_items.load() == 0; };

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

        try {
            m_array.reset(new std::atomic<tKeyType>[m_arraySize]()); //value init
        } catch (std::bad_alloc &e) {
            std::cerr << e.what() << std::endl;
            raise(SIGINT);
        }

        tbb::parallel_for(std::size_t(0), oldSize, [&oldArray, this](const uint i) {
            if (oldArray[i].load(std::memory_order_relaxed) != 0)
                insert(oldArray[i].load(std::memory_order_relaxed));
        });
    }

    void unsafe_merge(Concurrent_LP_Set &&other) {

        unsafe_rehash((capacity() + other.capacity()));

        tbb::parallel_for(std::size_t(0), other.capacity(), [&other, this](const uint i) {
            if (other.m_array[i].load(std::memory_order_relaxed) != 0)
                insert(other.m_array[i].load(std::memory_order_relaxed));
        });
    }


    template<class Set>
    void unsafe_merge(Concurrent_LP_Set &&other, const Set &filter) {

        VTUNE_TASK(MergeAllocate);
        std::size_t oldSize = m_arraySize;
        m_arraySize = nextPow2((capacity() + other.capacity()));
        m_hasher.l = log2(m_arraySize);

        auto oldArray = std::move(m_array);
        m_items.store(0, std::memory_order_relaxed);

        try {
            m_array.reset(new std::atomic<tKeyType>[m_arraySize]); //random init
        } catch (std::bad_alloc &e) {
            std::cerr << e.what() << std::endl;
            raise(SIGINT);
        }
        VTUNE_END_TASK(MergeAllocate);

        VTUNE_TASK(MergeFill);
        tbb::parallel_for(tbb::blocked_range<std::size_t>(0, m_arraySize), [this](const auto &r) {
            for (auto i = r.begin(); i != r.end(); ++i) {
                m_array[i].store(0, std::memory_order_relaxed);
            }
        });
        VTUNE_END_TASK(MergeFill);

        VTUNE_TASK(MergeCopy);
        tbb::parallel_for(tbb::blocked_range<std::size_t>(0, std::max(oldSize, other.capacity())),
                          [&oldArray, oldSize, &other, &filter, this](const auto &r) {

                              tKeyType val = 0;
                              for (auto i = r.begin(); i != r.end(); ++i) {
                                  if (i < oldSize) {
                                      val = oldArray[i].load(std::memory_order_relaxed);

                                      if (val != 0 && !filter.count(val))
                                          this->insert(val);
                                  }

                                  if (i < other.capacity()) {
                                      val = other.m_array[i].load(std::memory_order_relaxed);

                                      if (val != 0 && !filter.count(val))
                                          this->insert(val);
                                  }
                              }
                          });
    }

private:
    void helpGrowing() const {
        std::unique_lock<std::mutex> lk(m_growWaitMtx);
        m_growCV.wait(lk, [this] { return !m_fGrowing.load(); });
    }

    void grow() {

        bool exp = false;
        if (m_fGrowing.compare_exchange_strong(exp, true)) {

            VTUNE_TASK(GrowAllocate);
            std::size_t newSize = nextPow2(m_arraySize << 1);

            auto newHasher(m_hasher);
            newHasher.l = log2(newSize);

            std::unique_ptr<std::atomic<tKeyType>[]> newArray;
            try {
                newArray.reset(new std::atomic<tKeyType>[newSize]); //random init
            } catch (std::bad_alloc &e) {
                std::cerr << e.what() << std::endl;
                raise(SIGINT);
            }
            VTUNE_END_TASK(GrowAllocate);

            VTUNE_TASK(GrowFill);
            tbb::parallel_for(tbb::blocked_range<std::size_t>(0, newSize), [&newArray](const auto &r) {
                for (auto i = r.begin(); i != r.end(); ++i) {
                    newArray[i].store(0, std::memory_order_relaxed);
                }
            });
            VTUNE_END_TASK(GrowFill);

            VTUNE_TASK(GrowCopy);
            tbb::parallel_for(std::size_t(0), m_arraySize, [&newArray, &newSize, &newHasher, this](const uint i) {
                auto val = m_array[i].load(std::memory_order_relaxed);
                if (val != 0)
                    migrate(val, newArray, newSize, newHasher);
            });
            VTUNE_END_TASK(GrowCopy);

            //TODO this is all one transaction
            m_array = std::move(newArray);
            m_arraySize = newSize;
            m_hasher.l = newHasher.l;

            //wake up other threads;
            m_fGrowing.store(false);
            m_growCV.notify_all();
        } else {
            //somebody else already locked the mutex
            helpGrowing();
        }
    }

    void migrate(const tKeyType &key, std::unique_ptr<std::atomic<tKeyType>[]> &array,
                 const std::size_t &size, const Hasher<tKeyType> &hasher) {

        ASSERT(key != 0);

        for (tKeyType idx = hasher(key); ; idx++) {
            idx &= size - 1;
            ASSERT(idx < size);

            // Load the key that was there.
            tKeyType probedKey = array[idx].load();

            if (probedKey == key) {
                return; // the key is already in the set, return false;
            }
            else {
                // The entry was either free, or contains another key.
                if (probedKey != 0)
                    continue; // Usually, it contains another key. Keep probing.
                // The entry was free. Now let's try to take it using a CAS.
                tKeyType prevKey = 0;
                bool cas = array[idx].compare_exchange_strong(prevKey, key);
                if (cas) {
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

    std::atomic<bool> m_fGrowing;
    mutable std::mutex m_growWaitMtx;
    mutable std::condition_variable m_growCV;

};
