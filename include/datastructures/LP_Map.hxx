#pragma once

#include "LP_Set.hxx"

class LP_Map {

public:
    typedef uint tKeyType;
    typedef uint tValueType;

public:

    LP_Map(std::size_t size, Hasher<tKeyType> hasher = Hasher<tKeyType>())
            : m_items(0),
              m_keys(nullptr),
              m_hasher(hasher) {
        // Initialize cells
        m_arraySize = nextPow2(size);

        m_keys = std::make_unique<std::vector<tKeyType>>();
        m_keys->resize(m_arraySize);

        m_values = std::make_unique<std::vector<tValueType>>();
        m_values->resize(m_arraySize);


        m_hasher.l = log2(m_arraySize);
    }

    LP_Map(LP_Map &&other)
            : m_arraySize(other.m_arraySize),
              m_items(other.m_items),
              m_keys(std::move(other.m_keys)),
              m_values(std::move(other.m_values)),
              m_hasher(std::move(other.m_hasher)) { }

    LP_Map &operator=(LP_Map &&other) {
        m_arraySize = other.m_arraySize;
        m_items = other.m_items;
        m_keys = std::move(other.m_keys);
        m_values = std::move(other.m_values);
        m_hasher = std::move(other.m_hasher);

        return *this;
    }

    ~LP_Map() { }

    bool insert(const tKeyType &key, const tValueType &value) {
        ASSERT(key != 0);

        if (!m_rehashing && m_items / m_arraySize > 0.5)
            rehash(m_arraySize << 1);

        for (tKeyType idx = m_hasher(key); ; idx++) {
            idx &= m_arraySize - 1;
            ASSERT(idx < m_arraySize);

            // Load the key that was there.
            tKeyType probedKey = m_keys->at(idx);

            if (probedKey == key) {
                m_values->at(idx) = value; // key already exists, update value
                return false; // we didn't insert a new key
            }
            else {
                // The entry was either free, or contains another key.
                if (probedKey != 0)
                    continue; // Usually, it contains another key. Keep probing.
                // The entry was free. take it
                m_keys->at(idx) = key;
                m_values->at(idx) = value;
                ++m_items;
                return true;
            }
        }
    }

    bool insert(const std::pair<tKeyType, tValueType> & pair) {
        return insert(pair.first, pair.second);
    }

    tValueType get(const tKeyType &key) const {
        ASSERT(key != 0);

        for (tKeyType idx = m_hasher(key); ; idx++) {
            idx &= m_arraySize - 1;
            tKeyType probedKey = m_keys->at(idx);
            if (probedKey == key)
                return m_values->at(idx);;
            if (probedKey == 0)
                return 0;
        }

    }

    bool contains(const tKeyType &key) const {
        return get(key) != 0;
    }

    std::size_t count(const tKeyType &key) const {
        return contains(key);
    }


    bool empty() const { return m_items == 0; };

    bool empty(const std::size_t idx) const {
        return m_keys->at(idx) == 0;
    };

    std::pair<tKeyType, tValueType> at(const std::size_t idx) const {

        return std::make_pair(m_keys->at(idx), m_values->at(idx));
    }

    auto capacity() const { return m_arraySize; }

    auto size() const { return m_items; }

    void rehash(std::size_t newSize) {

        m_rehashing = true;

        std::size_t oldSize = m_arraySize;
        m_arraySize = nextPow2(newSize);
        m_hasher.l = log2(m_arraySize);

        auto oldKeys = std::move(m_keys);
        auto oldValues = std::move(m_values);
        m_items = 0;

        m_keys = std::make_unique<std::vector<tKeyType>>();
        m_keys->resize(m_arraySize);

        m_values = std::make_unique<std::vector<tValueType>>();
        m_values->resize(m_arraySize);

        for (std::size_t i = 0; i < oldSize; ++i) {
            if (oldKeys->at(i) != 0)
                insert(oldKeys->at(i), oldValues->at(i));
        }

        m_rehashing = false;

    }

    void merge(LP_Map &&other) {

        rehash((capacity() + other.capacity()) << 1);

        for (std::size_t i = 0; i < other.capacity(); ++i) {
            if (other.m_keys->at(i) != 0)
                insert(other.m_keys->at(i), other.m_values->at(i));
        }
    }

    template<class Set>
    void rehash(std::size_t newSize, const Set &filter) {

        m_rehashing = true;

        std::size_t oldSize = m_arraySize;
        m_arraySize = nextPow2(newSize);
        m_hasher.l = log2(m_arraySize);

        auto oldKeys = std::move(m_keys);
        auto oldValues = std::move(m_values);
        m_items = 0;

        m_keys = std::make_unique<std::vector<tKeyType>>();
        m_keys->resize(m_arraySize);

        m_values = std::make_unique<std::vector<tValueType>>();
        m_values->resize(m_arraySize);

        for (std::size_t i = 0; i < oldSize; ++i) {
            if (oldKeys->at(i) != 0 && !filter.count(oldKeys->at(i)))
                insert(oldKeys->at(i), oldValues->at(i));
        }

        m_rehashing = false;
    }

    template<class Set>
    void merge(LP_Map &&other, const Set &filter) {
        rehash((capacity() + other.capacity()) << 1, filter);

        for (std::size_t i = 0; i < other.capacity(); ++i) {
            if (other.m_keys->at(i) != 0 && !filter.count(other.m_keys->at(i)))
                insert(other.m_keys->at(i), other.m_values->at(i));
        }
    }


public:
    typedef _detail::iterator<LP_Map, std::pair<const tKeyType, tValueType>> iterator;
    typedef _detail::range_type<LP_Map, iterator> range_type;

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
    std::unique_ptr<std::vector<tKeyType>> m_keys;
    std::unique_ptr<std::vector<tValueType>> m_values;
    Hasher<tKeyType> m_hasher;

    bool m_rehashing = false;

};

class Concurrent_LP_Map {

public:
    typedef uint tKeyType;
    typedef uint tValueType;

public:

    Concurrent_LP_Map(std::size_t size, Hasher<tKeyType> hasher = Hasher<tKeyType>())
            : m_items(0),
              m_keys(nullptr),
              m_hasher(hasher) {
        // Initialize cells
        m_arraySize = nextPow2(size);
        m_keys = std::unique_ptr<std::atomic<tKeyType>[]>(new std::atomic<tKeyType>[m_arraySize]()); // zero init
        m_values = std::unique_ptr<std::atomic<tKeyType>[]>(new std::atomic<tValueType>[m_arraySize]); // random init

        m_hasher.l = log2(m_arraySize);
    }

    Concurrent_LP_Map(Concurrent_LP_Map &&other)
            : m_arraySize(other.m_arraySize),
              m_items(other.m_items.load(std::memory_order_relaxed)),
              m_keys(std::move(other.m_keys)),
              m_values(std::move(other.m_values)),
              m_hasher(std::move(other.m_hasher)) { }

    Concurrent_LP_Map &operator=(Concurrent_LP_Map &&other) {
        m_arraySize = other.m_arraySize;
        m_items.store(other.m_items.load(std::memory_order_relaxed), std::memory_order_relaxed);
        m_keys = std::move(other.m_keys);
        m_values = std::move(other.m_values);
        m_hasher = std::move(other.m_hasher);

        return *this;
    }

    ~Concurrent_LP_Map() { }

    bool insert(const tKeyType &key, const tValueType &value) {
        ASSERT(key != 0);

        if (m_items.load(std::memory_order_relaxed) / m_arraySize > 0.75)
            throw std::length_error("Overfull Concurrent Map");

        for (tKeyType idx = m_hasher(key); ; idx++) {
            idx &= m_arraySize - 1;
            ASSERT(idx < m_arraySize);

            // Load the key that was there.
            tKeyType probedKey = m_keys[idx].load(std::memory_order_relaxed);

            if (probedKey == key) {
                m_values[idx].store(value, std::memory_order_relaxed); // update value
                return false; // the key is already in the set, return false;
            }
            else {
                // The entry was either free, or contains another key.
                if (probedKey != 0)
                    continue; // Usually, it contains another key. Keep probing.
                // The entry was free. Now let's try to take it using a CAS.
                tKeyType prevKey = 0;
                bool cas = m_keys[idx].compare_exchange_strong(prevKey, key, std::memory_order_relaxed);
                if (cas) {
                    m_values[idx].store(value, std::memory_order_relaxed); // insert value
                    ++m_items;
                    return true; // we just added the key to the set
                }
                else if (prevKey == key) {
                    m_values[idx].store(value, std::memory_order_relaxed); // update value
                    return false; // the key was already added by another thread update it
                } else
                    continue; // another thread inserted a different key in this position
            }
        }
    }

    bool insert(const std::pair<tKeyType, tValueType> & pair) {
        return insert(pair.first, pair.second);
    }

    tValueType get(const tKeyType &key) const {
        ASSERT(key != 0);
        for (tKeyType idx = m_hasher(key); ; idx++) {
            idx &= m_arraySize - 1;
            tKeyType probedKey = m_keys[idx].load(std::memory_order_relaxed);
            if (probedKey == key)
                return m_values[idx].load(std::memory_order_relaxed);;
            if (probedKey == 0)
                return 0;
        }

    }

    bool contains(const tKeyType &key) const {
        return get(key) != 0;
    }

    std::size_t count(const tKeyType &key) const {
        return contains(key);
    }

    bool empty() const { return m_items.load(std::memory_order_relaxed) == 0; };

    bool empty(const std::size_t idx) const { return m_keys[idx].load(std::memory_order_relaxed) == 0; };

    std::pair<tKeyType, tValueType> at(const std::size_t idx) const {
        return std::make_pair(m_keys[idx].load(std::memory_order_relaxed),
                              m_values[idx].load(std::memory_order_relaxed));
    }

    auto capacity() const { return m_arraySize; }

    auto size() const { return m_items.load(); }

    void unsafe_rehash(std::size_t newSize) {
        std::size_t oldSize = m_arraySize;
        m_arraySize = nextPow2(newSize);
        m_hasher.l = log2(m_arraySize);

        auto oldKeys = std::move(m_keys);
        auto oldValues = std::move(m_values);
        m_items.store(0, std::memory_order_relaxed);

        m_keys = std::unique_ptr<std::atomic<tKeyType>[]>(new std::atomic<tKeyType>[m_arraySize]()); //zero init
        m_values = std::unique_ptr<std::atomic<tValueType>[]>(new std::atomic<tKeyType>[m_arraySize]); //random init

        tbb::parallel_for(std::size_t(0), oldSize, [&oldKeys, &oldValues, this](const uint i) {
            if (oldKeys[i].load(std::memory_order_relaxed) != 0)
                insert(oldKeys[i].load(std::memory_order_relaxed), oldValues[i].load(std::memory_order_relaxed));
        });
    }

    void unsafe_merge(Concurrent_LP_Map &&other) {

        unsafe_rehash((capacity() + other.capacity()) << 1);

        tbb::parallel_for(std::size_t(0), other.capacity(), [&other, this](const uint i) {
            if (other.m_keys[i].load(std::memory_order_relaxed) != 0)
                insert(other.m_keys[i].load(std::memory_order_relaxed),
                       other.m_values[i].load(std::memory_order_relaxed));
        });
    }


    template<class Set>
    void unsafe_merge(Concurrent_LP_Map &&other, const Set &filter) {

        std::size_t oldSize = m_arraySize;
        m_arraySize = nextPow2((capacity() + other.capacity()) << 1);
        m_hasher.l = log2(m_arraySize);

        auto oldKeys = std::move(m_keys);
        auto oldValues = std::move(m_values);
        m_items.store(0, std::memory_order_relaxed);

        m_keys = std::unique_ptr<std::atomic<tKeyType>[]>(new std::atomic<tKeyType>[m_arraySize]()); // zero init
        m_keys = std::unique_ptr<std::atomic<tValueType>[]>(new std::atomic<tKeyType>[m_arraySize]); // random init

        tbb::parallel_for(tbb::blocked_range<std::size_t>(0, std::max(oldSize, other.capacity())),
                          [&oldKeys, &oldValues, oldSize, &other, &filter, this](const auto &r) {

                              tKeyType val = 0;
                              for (auto i = r.begin(); i != r.end(); ++i) {
                                  if (i < oldSize) {
                                      val = oldKeys[i].load(std::memory_order_relaxed);

                                      if (val != 0 && !filter.count(val))
                                          this->insert(val, oldValues[i].load(std::memory_order_relaxed));
                                  }

                                  if (i < other.capacity()) {
                                      val = other.m_keys[i].load(std::memory_order_relaxed);

                                      if (val != 0 && !filter.count(val))
                                          this->insert(val, other.m_values[i].load(std::memory_order_relaxed));
                                  }
                              }
                          });
    }

public:
    typedef _detail::iterator<Concurrent_LP_Map, std::pair<const tKeyType, tValueType>> iterator;
    typedef _detail::range_type<Concurrent_LP_Map, iterator> range_type;

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
    std::unique_ptr<std::atomic<tKeyType>[]> m_keys;
    std::unique_ptr<std::atomic<tValueType>[]> m_values;
    Hasher<tKeyType> m_hasher;

};
