#include <atomic>
#include <functional>
#include <cstring>

#include <iostream>
#include <iomanip>

#include "Misc.h"
#include "ASSERT.h"
#include "Random.h"

template<typename T>
struct Hasher {

    Hasher() {
        std::uniform_int_distribution<T> dist(std::numeric_limits<T>::min(), std::numeric_limits<T>::max());
        while (a = dist(startGen), a % 2 == 0) { } // find odd value

        l = sizeof(T) * CHAR_BIT;
    }

    T operator()(const T x) {
        return (a * x) >> (sizeof(T) * CHAR_BIT - l);
    }

    T a;
    T l;
};

class LP_Set {

public:
    typedef uint tKeyType;

public:

    LP_Set(std::size_t size, Hasher<tKeyType> &hasher)
            : m_items(0),
              m_array(nullptr),
              m_hasher(hasher) {
        // Initialize cells
        m_arraySize = nextPow2(size);
        m_array = new tKeyType[m_arraySize];
        std::fill(m_array, m_array + m_arraySize, 0);

        m_hasher.l = log2(m_arraySize);
    }

    ~LP_Set() {
        if (m_array)
            delete[] m_array;
    }

    bool insert(const tKeyType &key) {
        ASSERT(key != 0);

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

    bool contains(const tKeyType &key) {
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

    bool erase(const tKeyType &key) {
        ASSERT(key != 0);
        for (tKeyType idx = m_hasher(key); ; idx++) {
            idx &= m_arraySize - 1;
            tKeyType probedKey = m_array[idx];
            if (probedKey == key) {
                m_array[idx] = 0;
                --m_items;
                return true; // we just deleted the key from the set
            } else if (probedKey == 0)
                return false;
        }

    }

    auto capacity() const { return m_arraySize; }

    auto size() const { return m_items; }

    void rehash(std::size_t newSize) {
        std::size_t oldSize = m_arraySize;
        m_arraySize = nextPow2(newSize);

        tKeyType *oldArray = m_array;

        m_array = new tKeyType[m_arraySize];
        std::fill(m_array, m_array + m_arraySize, 0);

        for (std::size_t i = 0; i < oldSize; ++i) {
            insert(oldArray[i]);
        }

        delete[] oldArray;
    }

    void merge(LP_Set &&other) {

        rehash(capacity() + other.capacity());

        for (std::size_t i = 0; i < other.capacity(); ++i) {
            insert(other.m_array[i]);
        }
    }


private:
    std::size_t m_arraySize;
    std::size_t m_items;
    tKeyType *m_array;
    Hasher<tKeyType> m_hasher;

};

class Concurrent_LP_Set {

public:
    typedef uint tKeyType;

public:

    Concurrent_LP_Set(std::size_t size, Hasher<tKeyType> &hasher)
            : m_items(0),
              m_array(nullptr),
              m_hasher(hasher) {
        // Initialize cells
        m_arraySize = nextPow2(size);
        m_array = new std::atomic<tKeyType>[m_arraySize];
        std::fill(m_array, m_array + m_arraySize, 0);

        m_hasher.l = log2(m_arraySize);
    }

    ~Concurrent_LP_Set() {
        if (m_array)
            delete[] m_array;
    }

    bool insert(const tKeyType &key) {
        ASSERT(key != 0);

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

    bool contains(const tKeyType &key) {
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

    bool erase(const tKeyType &key) {
        ASSERT(key != 0);
        for (tKeyType idx = m_hasher(key); ; idx++) {
            idx &= m_arraySize - 1;
            tKeyType probedKey = m_array[idx].load(std::memory_order_relaxed);
            if (probedKey == key) {
                tKeyType prevKey = key; // we found the key
                bool cas = m_array[idx].compare_exchange_strong(prevKey, 0, std::memory_order_relaxed);
                if (cas) {
                    --m_items;
                    return true; // we just deleted the key from the set
                } else
                    return false; //another thread must have deleted it
            } else if (probedKey == 0)
                return false;
        }

    }

    auto capacity() const { return m_arraySize; }

    auto size() const { return m_items.load(); }

    void unsafe_rehash(std::size_t newSize) {
        std::size_t oldSize = m_arraySize;
        m_arraySize = nextPow2(newSize);

        std::atomic<tKeyType> *oldArray = m_array;

        m_array = new std::atomic<tKeyType>[m_arraySize];
        std::fill(m_array, m_array + m_arraySize, 0);

        for (std::size_t i = 0; i < oldSize; ++i) {
            insert(oldArray[i].load(std::memory_order_relaxed));
        }

        delete[] oldArray;
    }

    void unsafe_merge(Concurrent_LP_Set &&other) {

        unsafe_rehash(capacity() + other.capacity());

        for (std::size_t i = 0; i < other.capacity(); ++i) {
            insert(other.m_array[i].load(std::memory_order_relaxed));
        }
    }


private:
    std::size_t m_arraySize;
    std::atomic<std::size_t> m_items;
    std::atomic<tKeyType> *m_array;
    Hasher<tKeyType> m_hasher;

};
