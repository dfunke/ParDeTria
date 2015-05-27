#pragma once

#include "LP_Set.hxx"

namespace _detail {

    template<class Container, typename KeyType, typename ValueType>
    struct filtered_iterator
            : public std::iterator<std::forward_iterator_tag, std::pair<const KeyType, ValueType>> {

    public:
        filtered_iterator(const Container &container, std::size_t idx, std::size_t bound, const KeyType &filter)
                : m_container(container),
                  m_idx(idx),
                  m_bound(bound),
                  m_filter(filter) { }

        filtered_iterator(const filtered_iterator &o)
                : m_container(o.m_container),
                  m_idx(o.m_idx),
                  m_bound(o.m_bound),
                  m_filter(o.m_filter) { }

        filtered_iterator &operator++() {
            while (m_idx = (m_idx + 1) & (m_container.capacity() - 1),
                    (m_idx != m_bound &&
                     (m_container.empty(m_idx) || m_container.at(m_idx).first != m_filter))) { }
            return *this;
        }

        filtered_iterator operator++(int) {
            auto tmp = *this;
            ++(*this);
            return tmp;
        }

        std::size_t operator-(const filtered_iterator &other) const {
            return m_idx - other.m_idx;
        }

        std::size_t operator+(const filtered_iterator &other) const {
            return m_idx + other.m_idx;
        }

        bool operator==(const filtered_iterator &j) const {
            return m_idx == j.m_idx;
        }

        bool operator!=(const filtered_iterator &j) const { return !(*this == j); }

        bool operator<(const filtered_iterator &other) const {
            return m_idx < other.m_idx;
        }

        bool operator<=(const filtered_iterator &other) const {
            return m_idx <= other.m_idx;
        }

        bool operator>(const filtered_iterator &other) const {
            return m_idx > other.m_idx;
        }

        bool operator>=(const filtered_iterator &other) const {
            return m_idx >= other.m_idx;
        }

        const auto operator*() const { return m_container.at(m_idx); }

        //const auto *operator->() const { return &m_container.at(m_idx); }

    protected:
        const Container &m_container;
        std::size_t m_idx;
        std::size_t m_bound;
        const KeyType m_filter;
    };
}

class LP_MultiMap {

public:
    typedef uint tKeyType;
    typedef uint tValueType;

public:
    typedef _detail::iterator<LP_MultiMap, std::pair<const tKeyType, tValueType>> iterator;
    typedef _detail::filtered_iterator<LP_MultiMap, tKeyType, tValueType> filtered_iterator;
    typedef _detail::range_type<LP_MultiMap, iterator> range_type;

public:

    LP_MultiMap(std::size_t size, Hasher<tKeyType> hasher = Hasher<tKeyType>())
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

    LP_MultiMap(LP_MultiMap &&other)
            : m_arraySize(other.m_arraySize),
              m_items(other.m_items),
              m_keys(std::move(other.m_keys)),
              m_values(std::move(other.m_values)),
              m_hasher(std::move(other.m_hasher)) { }

    LP_MultiMap &operator=(LP_MultiMap &&other) {
        m_arraySize = other.m_arraySize;
        m_items = other.m_items;
        m_keys = std::move(other.m_keys);
        m_values = std::move(other.m_values);
        m_hasher = std::move(other.m_hasher);

        return *this;
    }

    ~LP_MultiMap() { }

    bool insert(const tKeyType &key, const tValueType &value) {

        ASSERT(key != 0);

        if (!m_rehashing && m_items / m_arraySize > 0.5)
            rehash(m_arraySize << 1);

        for (tKeyType idx = m_hasher(key); ; idx++) {
            idx &= m_arraySize - 1;
            ASSERT(idx < m_arraySize);

            // Load the key that was there.
            tKeyType probedKey = m_keys->at(idx);

            if (probedKey != 0)
                continue; // Usually, it contains another key. Keep probing.

            // The entry was free. take it
            m_keys->at(idx) = key;
            m_values->at(idx) = value;
            ++m_items;
            return true;
        }
    }

    bool insert(const std::pair<tKeyType, tValueType> &pair) {
        return insert(pair.first, pair.second);
    }

    std::pair<filtered_iterator, filtered_iterator> get(const tKeyType &key) const {
        ASSERT(key != 0);

        std::size_t firstIdx = m_hasher(key) & (m_arraySize - 1);
        std::size_t lastIdx = firstIdx;

        for (tKeyType idx = firstIdx; ; idx++) {
            idx &= m_arraySize - 1;
            tKeyType probedKey = m_keys->at(idx);

            if (probedKey == 0 && idx == firstIdx)
                // the key can not be in the table, return to end() iterators
                return std::make_pair(filtered_iterator(*this, m_arraySize, m_arraySize, key),
                                      filtered_iterator(*this, m_arraySize, m_arraySize, key));
            else {
                if (probedKey != key && idx == firstIdx) {
                    //there are entries with the hash corresponding to our key, but they don't equal the key
                    //and we haven't found the first valid entry yet
                    ++firstIdx;
                    firstIdx &= m_arraySize - 1;
                }

                ++lastIdx;
                lastIdx &= m_arraySize - 1;

                if (probedKey == 0 /* implicitly && idx != firstIdx */)
                    // we have found the end of the chain of equal hashes
                    return std::make_pair(filtered_iterator(*this, firstIdx, lastIdx, key),
                                          filtered_iterator(*this, lastIdx, lastIdx, key));
            }
        }

    }

    bool contains(const tKeyType &key) const {
        ASSERT(key != 0);

        for (tKeyType idx = m_hasher(key); ; idx++) {
            idx &= m_arraySize - 1;
            tKeyType probedKey = m_keys->at(idx);
            if (probedKey == key)
                return true;
            if (probedKey == 0)
                return false;
        }
    }

    std::size_t count(const tKeyType &key) const {
        return contains(key);
    }


    bool empty() const { return m_items == 0; };

    bool empty(const std::size_t idx) const {
        return m_keys->at(idx) == 0;
    };

    std::pair<const tKeyType, tValueType> at(const std::size_t idx) const {

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

    void merge(LP_MultiMap &&other) {

        rehash((capacity() + other.capacity()));

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
    void merge(LP_MultiMap &&other, const Set &filter) {
        rehash((capacity() + other.capacity()), filter);

        for (std::size_t i = 0; i < other.capacity(); ++i) {
            if (other.m_keys->at(i) != 0 && !filter.count(other.m_keys->at(i)))
                insert(other.m_keys->at(i), other.m_values->at(i));
        }
    }

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

class Concurrent_LP_MultiMap {

public:
    typedef uint tKeyType;
    typedef uint tValueType;

public:
    typedef _detail::iterator<Concurrent_LP_MultiMap, std::pair<const tKeyType, tValueType>> iterator;
    typedef _detail::filtered_iterator<Concurrent_LP_MultiMap, tKeyType, tValueType> filtered_iterator;
    typedef _detail::range_type<Concurrent_LP_MultiMap, iterator> range_type;

public:

    Concurrent_LP_MultiMap(std::size_t size, Hasher<tKeyType> hasher = Hasher<tKeyType>())
            : m_items(0),
              m_keys(nullptr),
              m_hasher(hasher) {
        // Initialize cells
        m_arraySize = nextPow2(size);

        try {
            m_keys = std::unique_ptr<std::atomic<tKeyType>[]>(new std::atomic<tKeyType>[m_arraySize]()); // zero init
            m_values = std::unique_ptr<std::atomic<tKeyType>[]>(
                    new std::atomic<tValueType>[m_arraySize]); // random init
        } catch (std::bad_alloc &e) {
            std::cerr << e.what() << std::endl;
            raise(SIGINT);
        }

        m_hasher.l = log2(m_arraySize);
    }

    Concurrent_LP_MultiMap(Concurrent_LP_MultiMap &&other)
            : m_arraySize(other.m_arraySize),
              m_items(other.m_items.load()),
              m_keys(std::move(other.m_keys)),
              m_values(std::move(other.m_values)),
              m_hasher(std::move(other.m_hasher)) { }

    Concurrent_LP_MultiMap &operator=(Concurrent_LP_MultiMap &&other) {
        m_arraySize = other.m_arraySize;
        m_items.store(other.m_items.load());
        m_keys = std::move(other.m_keys);
        m_values = std::move(other.m_values);
        m_hasher = std::move(other.m_hasher);

        return *this;
    }

    ~Concurrent_LP_MultiMap() { }

    bool insert(const tKeyType &key, const tValueType &value) {
        ASSERT(key != 0);


        std::size_t steps = 0; //count number of steps
        for (tKeyType idx = m_hasher(key); ; idx++) {

            if (steps > m_arraySize / 4) {
                throw std::length_error("Overfull MultiMap");
                /*//we stepped through more than a quarter the array > grow
                grow();

                //reset insertion
                steps = 0;
                idx = m_hasher(key);*/
            }

            idx &= m_arraySize - 1;
            ASSERT(idx < m_arraySize);

            // Load the key that was there.
            tKeyType probedKey = m_keys[idx].load();

            if (probedKey != 0) {
                ++steps;
                continue; // Usually, it contains another key. Keep probing.
            }
            // The entry was free. Now let's try to take it using a CAS.
            tKeyType prevKey = 0;
            bool cas = m_keys[idx].compare_exchange_strong(prevKey, key);
            if (cas) {
                m_values[idx].store(value, std::memory_order_relaxed); // insert value
                ++m_items;
                return true; // we just added the key to the set
            } else {
                ++steps;
                continue; // another thread inserted in this position
            }
        }
    }

    bool insert(const std::pair<tKeyType, tValueType> &pair) {
        return insert(pair.first, pair.second);
    }

    std::pair<filtered_iterator, filtered_iterator> get(const tKeyType &key) const {
        ASSERT(key != 0);

        std::size_t firstIdx = m_hasher(key) & (m_arraySize - 1);
        std::size_t lastIdx = firstIdx;

        for (tKeyType idx = firstIdx; ; idx++) {
            idx &= m_arraySize - 1;
            tKeyType probedKey = m_keys[idx].load();

            if (probedKey == 0 && idx == firstIdx)
                // the key can not be in the table, return to end() iterators
                return std::make_pair(filtered_iterator(*this, m_arraySize, m_arraySize, key),
                                      filtered_iterator(*this, m_arraySize, m_arraySize, key));
            else {
                if (probedKey != key && idx == firstIdx) {
                    //there are entries with the hash corresponding to our key, but they don't equal the key
                    //and we haven't found the first valid entry yet
                    ++firstIdx;
                    firstIdx &= m_arraySize - 1;
                }

                ++lastIdx;
                lastIdx &= m_arraySize - 1;

                if (probedKey == 0 /* implicitly && idx != firstIdx */)
                    // we have found the end of the chain of equal hashes
                    return std::make_pair(filtered_iterator(*this, firstIdx, lastIdx, key),
                                          filtered_iterator(*this, lastIdx, lastIdx, key));
            }
        }

    }

    bool contains(const tKeyType &key) const {
        ASSERT(key != 0);
        for (tKeyType idx = m_hasher(key); ; idx++) {
            idx &= m_arraySize - 1;
            tKeyType probedKey = m_keys[idx].load();
            if (probedKey == key)
                return true;
            if (probedKey == 0)
                return false;
        }

    }

    std::size_t count(const tKeyType &key) const {
        return contains(key);
    }

    bool empty() const { return m_items.load() == 0; };

    bool empty(const std::size_t idx) const { return m_keys[idx].load(std::memory_order_relaxed) == 0; };

    std::pair<const tKeyType, tValueType> at(const std::size_t idx) const {
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

        try {
            m_keys = std::unique_ptr<std::atomic<tKeyType>[]>(new std::atomic<tKeyType>[m_arraySize]()); //zero init
            m_values = std::unique_ptr<std::atomic<tValueType>[]>(new std::atomic<tKeyType>[m_arraySize]); //random init
        } catch (std::bad_alloc &e) {
            std::cerr << e.what() << std::endl;
            raise(SIGINT);
        }

        tbb::parallel_for(std::size_t(0), oldSize, [&oldKeys, &oldValues, this](const uint i) {
            if (oldKeys[i].load(std::memory_order_relaxed) != 0)
                insert(oldKeys[i].load(std::memory_order_relaxed), oldValues[i].load(std::memory_order_relaxed));
        });
    }

    void unsafe_merge(Concurrent_LP_MultiMap &&other) {

        unsafe_rehash((capacity() + other.capacity()));

        tbb::parallel_for(std::size_t(0), other.capacity(), [&other, this](const uint i) {
            if (other.m_keys[i].load(std::memory_order_relaxed) != 0)
                insert(other.m_keys[i].load(std::memory_order_relaxed),
                       other.m_values[i].load(std::memory_order_relaxed));
        });
    }


    template<class Set>
    void unsafe_merge(Concurrent_LP_MultiMap &&other, const Set &filter) {

        std::size_t oldSize = m_arraySize;
        m_arraySize = nextPow2((capacity() + other.capacity()));
        m_hasher.l = log2(m_arraySize);

        auto oldKeys = std::move(m_keys);
        auto oldValues = std::move(m_values);
        m_items.store(0, std::memory_order_relaxed);

        try {
            m_keys = std::unique_ptr<std::atomic<tKeyType>[]>(new std::atomic<tKeyType>[m_arraySize]()); // zero init
            m_keys = std::unique_ptr<std::atomic<tValueType>[]>(new std::atomic<tKeyType>[m_arraySize]); // random init
        } catch (std::bad_alloc &e) {
            std::cerr << e.what() << std::endl;
            raise(SIGINT);
        }


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
