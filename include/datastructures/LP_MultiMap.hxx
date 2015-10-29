#pragma once

#include "LP_Set.hxx"

namespace _detail {

    template<class Container, typename KeyType, typename ValueType>
    struct filtered_iterator
            : public std::iterator<std::forward_iterator_tag, std::pair<const KeyType, ValueType>> {

    public:
        filtered_iterator(Container &container, std::size_t idx, std::size_t bound, const KeyType &filter)
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

        std::size_t operator+(const std::size_t &a) const {
            return m_idx + a;
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

        auto operator*() const { return m_container.at(m_idx); }

        //const auto *operator->() const { return &m_container.at(m_idx); }

        void half(filtered_iterator &begin, filtered_iterator & end){
            _setIdx(begin + ((end - begin) / 2), true);
        }

    protected:
        Container &m_container;
        std::size_t m_idx;
        std::size_t m_bound;
        const KeyType m_filter;
    };
}

template<typename K, typename V>
class Concurrent_LP_MultiMap;

template<typename K, typename V>
class LP_MultiMap {

    friend class Concurrent_LP_MultiMap<K, V>;

public:
    typedef _detail::iterator<const LP_MultiMap<K, V>, std::pair<const K, V>> const_iterator;
    typedef _detail::filtered_iterator<const LP_MultiMap<K, V>, K, V> const_filtered_iterator;
    typedef _detail::range_type<const LP_MultiMap<K, V>, const_iterator> const_range_type;

    LP_MultiMap(std::size_t size, Hasher<K> hasher = Hasher<K>())
            : m_items(0),
              m_keys(nullptr),
              m_hasher(hasher) {
        // Initialize cells
        m_arraySize = nextPow2(size);

        m_keys.reset(new K[m_arraySize]()); // zero init
        m_values.reset(new V[m_arraySize]); // random init
    }

    LP_MultiMap(LP_MultiMap<K, V> &&other)
            : m_arraySize(other.m_arraySize),
              m_items(other.m_items),
              m_keys(std::move(other.m_keys)),
              m_values(std::move(other.m_values)),
              m_hasher(std::move(other.m_hasher)) { }

    LP_MultiMap &operator=(LP_MultiMap<K, V> &&other) {
        m_arraySize = other.m_arraySize;
        m_items = other.m_items;
        m_keys = std::move(other.m_keys);
        m_values = std::move(other.m_values);
        m_hasher = std::move(other.m_hasher);

        return *this;
    }

    //conversion
    LP_MultiMap(Concurrent_LP_MultiMap<K, V> &&other);

    bool insert(const K &key, const V &value) {

        ASSERT(key != 0);

        if (!m_rehashing && m_items / m_arraySize >= 0.5)
            rehash(m_arraySize << 1);

        for (K idx = m_hasher(key); ; idx++) {
            idx &= m_arraySize - 1;
            ASSERT(idx < m_arraySize);

            // Load the key that was there.
            K probedKey = m_keys[idx];

            if (probedKey != 0)
                continue; // Usually, it contains another key. Keep probing.

            // The entry was free. take it
            m_keys[idx] = key;
            m_values[idx] = value;
            ++m_items;
            return true;
        }
    }

    bool insert(const std::pair<K, V> &pair) {
        return insert(pair.first, pair.second);
    }

    std::pair<const_filtered_iterator, const_filtered_iterator> get(const K &key) const {
        ASSERT(key != 0);

        std::size_t firstIdx = m_hasher(key) & (m_arraySize - 1);
        std::size_t lastIdx = firstIdx;

        for (K idx = firstIdx; ; idx++) {
            idx &= m_arraySize - 1;
            K probedKey = m_keys[idx];

            if (probedKey == 0 && idx == firstIdx)
                // the key can not be in the table, return to end() iterators
                return std::make_pair(const_filtered_iterator(*this, m_arraySize, m_arraySize, key),
                                      const_filtered_iterator(*this, m_arraySize, m_arraySize, key));
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
                    return std::make_pair(const_filtered_iterator(*this, firstIdx, lastIdx, key),
                                          const_filtered_iterator(*this, lastIdx, lastIdx, key));
            }
        }

    }

    bool contains(const K &key) const {
        ASSERT(key != 0);

        for (K idx = m_hasher(key); ; idx++) {
            idx &= m_arraySize - 1;
            K probedKey = m_keys[idx];
            if (probedKey == key)
                return true;
            if (probedKey == 0)
                return false;
        }
    }

    std::size_t count(const K &key) const {
        return contains(key);
    }


    bool empty() const { return m_items == 0; };

    bool empty(const std::size_t idx) const {
        return m_keys[idx] == 0;
    };

    std::pair<const K, V> at(const std::size_t idx) const {

        return std::make_pair(m_keys[idx], m_values[idx]);
    }

    auto capacity() const { return m_arraySize; }

    auto size() const { return m_items; }

    void rehash(std::size_t newSize) {

        m_rehashing = true;

        std::size_t oldSize = m_arraySize;
        m_arraySize = nextPow2(newSize);

        auto oldKeys = std::move(m_keys);
        auto oldValues = std::move(m_values);
        m_items = 0;

        m_keys.reset(new K[m_arraySize]()); // zero init
        m_values.reset(new V[m_arraySize]); // random init

        for (std::size_t i = 0; i < oldSize; ++i) {
            if (oldKeys[i] != 0)
                insert(oldKeys[i], oldValues[i]);
        }

        m_rehashing = false;

    }

    void merge(LP_MultiMap<K, V> &&other) {

        rehash((capacity() + other.capacity()));

        for (std::size_t i = 0; i < other.capacity(); ++i) {
            if (other.m_keys[i] != 0)
                insert(other.m_keys[i], other.m_values[i]);
        }
    }

    template<class Set>
    void rehash(std::size_t newSize, const Set &filter) {

        m_rehashing = true;

        std::size_t oldSize = m_arraySize;
        m_arraySize = nextPow2(newSize);

        auto oldKeys = std::move(m_keys);
        auto oldValues = std::move(m_values);
        m_items = 0;

        m_keys.reset(new K[m_arraySize]()); // zero init
        m_values.reset(new V[m_arraySize]); // random init

        for (std::size_t i = 0; i < oldSize; ++i) {
            if (oldKeys[i] != 0 && !filter.count(oldKeys[i]))
                insert(oldKeys[i], oldValues[i]);
        }

        m_rehashing = false;
    }

    template<class Set>
    void merge(LP_MultiMap<K, V> &&other, const Set &filter) {
        rehash((capacity() + other.capacity()), filter);

        for (std::size_t i = 0; i < other.capacity(); ++i) {
            if (other.m_keys[i] != 0 && !filter.count(other.m_keys[i]))
                insert(other.m_keys[i], other.m_values[i]);
        }
    }

    const_iterator begin() const {
        return const_iterator(*this, 0, true);
    }

    const_iterator end() const {
        return const_iterator(*this, m_arraySize, false);
    }

    const_range_type range() const {
        return const_range_type(*this);
    }

private:
    friend class boost::serialization::access;

    template<class Archive>
    void save(Archive &ar, __attribute__((unused)) const unsigned int version) const {
        // invoke serialization of the base class
        ar << m_arraySize;
        ar << m_items;
        ar << boost::serialization::make_array<K>(m_keys.get(), m_arraySize);
        ar << boost::serialization::make_array<V>(m_values.get(), m_arraySize);
        ar << m_hasher.a;
    }

    template<class Archive>
    void load(Archive &ar, __attribute__((unused)) const unsigned int version) {
        // invoke serialization of the base class
        ar >> m_arraySize;
        ar >> m_items;

        try {
            m_keys.reset(new K[m_arraySize]); //random init
            m_values.reset(new V[m_arraySize]); //random init
        } catch (std::bad_alloc &e) {
            std::cerr << e.what() << std::endl;
            raise(SIGINT);
        }
        ar >> boost::serialization::make_array<K>(m_keys.get(), m_arraySize);
        ar >> boost::serialization::make_array<V>(m_values.get(), m_arraySize);
        ar >> m_hasher.a;
    }

    BOOST_SERIALIZATION_SPLIT_MEMBER()

private:
    std::size_t m_arraySize;
    std::size_t m_items;
    std::unique_ptr<K[]> m_keys;
    std::unique_ptr<V[]> m_values;
    Hasher<K> m_hasher;

    bool m_rehashing = false;

};

template<typename K, typename V>
class Concurrent_LP_MultiMap {

    friend class LP_MultiMap<K, V>;

    friend class GrowingHashTable<Concurrent_LP_MultiMap<K, V>>;

    friend class GrowingHashTableHandle<Concurrent_LP_MultiMap<K, V>>;

    friend class ConstGrowingHashTableHandle<Concurrent_LP_MultiMap<K, V>>;

public:
    typedef _detail::iterator<const Concurrent_LP_MultiMap<K, V>, std::pair<const K, V>> const_iterator;
    typedef _detail::filtered_iterator<const Concurrent_LP_MultiMap<K, V>, K, V> const_filtered_iterator;
    typedef _detail::range_type<const Concurrent_LP_MultiMap<K, V>, const_iterator> const_range_type;

public:

    Concurrent_LP_MultiMap(std::size_t size, const uint version = 0,
                           Hasher<K> hasher = Hasher<K>())
            : m_items(0),
              m_keys(nullptr),
              m_version(version),
              m_hasher(hasher),
              m_currentCopyBlock(0) {
        // Initialize cells
        m_arraySize = nextPow2(size);

        try {
            m_keys.reset(new std::atomic<K>[m_arraySize]()); // zero init
            m_values.reset(new std::atomic<V>[m_arraySize]); // random init
        } catch (std::bad_alloc &e) {
            std::cerr << e.what() << std::endl;
            raise(SIGINT);
        }

    }

    Concurrent_LP_MultiMap(Concurrent_LP_MultiMap<K, V> &&other)
            : m_arraySize(other.m_arraySize),
              m_items(other.m_items.load()),
              m_keys(std::move(other.m_keys)),
              m_values(std::move(other.m_values)),
              m_version(other.m_version),
              m_hasher(std::move(other.m_hasher)),
              m_currentCopyBlock(0) { }

    Concurrent_LP_MultiMap &operator=(Concurrent_LP_MultiMap<K, V> &&other) {
        m_arraySize = other.m_arraySize;
        m_items.store(other.m_items.load());
        m_keys = std::move(other.m_keys);
        m_values = std::move(other.m_values);
        m_version = other.m_version;
        m_hasher = std::move(other.m_hasher);
        m_currentCopyBlock.store(other.m_currentCopyBlock.load());

        return *this;
    }

    //conversion
    Concurrent_LP_MultiMap(LP_MultiMap<K, V> &&other);


    InsertReturn insert(const K &key, const V &value) {
        ASSERT(key != 0);


        if (m_items.load() >= m_arraySize >> 1)
            return InsertReturn::State::Full;

        std::size_t cols = 0; //count number of steps
        for (K idx = m_hasher(key); ; idx++) {

            //if (cols > c_growThreshold) {
            //    return InsertReturn::State::Full;
            //}

            idx &= m_arraySize - 1;
            ASSERT(idx < m_arraySize);

            // Load the key that was there.
            K probedKey = m_keys[idx].load();

            if (probedKey != 0) {
                cols += probedKey != key;
                continue; // Usually, it contains another key. Keep probing.
            }
            // The entry was free. Now let's try to take it using a CAS.
            K prevKey = 0;
            bool cas = m_keys[idx].compare_exchange_strong(prevKey, key);
            if (cas) {
                m_values[idx].store(value, std::memory_order_relaxed); // insert value
                ++m_items;
                return true; // we just added the key to the set
            } else {
                cols += prevKey != key;
                continue; // another thread inserted in this position
            }
        }
    }

    InsertReturn insert(const std::pair<K, V> &pair) {
        return insert(pair.first, pair.second);
    }

    std::pair<const_filtered_iterator, const_filtered_iterator> get(const K &key) const {
        ASSERT(key != 0);

        std::size_t firstIdx = m_hasher(key) & (m_arraySize - 1);
        std::size_t lastIdx = firstIdx;

        for (K idx = firstIdx; ; idx++) {
            idx &= m_arraySize - 1;
            K probedKey = m_keys[idx].load();

            if (probedKey == 0 && idx == firstIdx)
                // the key can not be in the table, return to end() iterators
                return std::make_pair(const_filtered_iterator(*this, m_arraySize, m_arraySize, key),
                                      const_filtered_iterator(*this, m_arraySize, m_arraySize, key));
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
                    return std::make_pair(const_filtered_iterator(*this, firstIdx, lastIdx, key),
                                          const_filtered_iterator(*this, lastIdx, lastIdx, key));
            }
        }

    }

    bool contains(const K &key) const {
        ASSERT(key != 0);
        for (K idx = m_hasher(key); ; idx++) {
            idx &= m_arraySize - 1;
            K probedKey = m_keys[idx].load();
            if (probedKey == key)
                return true;
            if (probedKey == 0)
                return false;
        }

    }

    std::size_t count(const K &key) const {
        return contains(key);
    }

    bool empty() const { return m_items.load() == 0; };

    bool empty(const std::size_t idx) const { return m_keys[idx].load(std::memory_order_relaxed) == 0; };

    std::pair<const K, V> at(const std::size_t idx) const {
        return std::make_pair(m_keys[idx].load(std::memory_order_relaxed),
                              m_values[idx].load(std::memory_order_relaxed));
    }

    auto capacity() const { return m_arraySize; }

    auto size() const { return m_items.load(); }

    void unsafe_rehash(std::size_t newSize) {
        std::size_t oldSize = m_arraySize;
        m_arraySize = nextPow2(newSize);

        auto oldKeys = std::move(m_keys);
        auto oldValues = std::move(m_values);
        m_items.store(0, std::memory_order_relaxed);

        try {
            m_keys.reset(new std::atomic<K>[m_arraySize]()); //zero init
            m_values.reset(new std::atomic<V>[m_arraySize]); //random init
        } catch (std::bad_alloc &e) {
            std::cerr << e.what() << std::endl;
            raise(SIGINT);
        }

        tbb::parallel_for(std::size_t(0), oldSize, [&oldKeys, &oldValues, this](const K i) {
            if (oldKeys[i].load(std::memory_order_relaxed) != 0)
                insert(oldKeys[i].load(std::memory_order_relaxed), oldValues[i].load(std::memory_order_relaxed));
        });
    }

    template<class Source>
    void unsafe_copy(const Source &other) {
        ASSERT(m_arraySize == other.size());
        for (std::size_t i = 0; i < m_arraySize; ++i) {
            auto pair = other.at(i);
            m_keys[i] = pair.first;
            m_values[i] = pair.second;
        }
    }

    template<class Source, class Filter>
    void unsafe_merge(Source &&other, const Filter &filter) {

        VTUNE_TASK(MultiMapMergeAllocate);
        std::size_t oldSize = m_arraySize;

        std::size_t combinedItems = size() + other.size();
        std::size_t newSize = nextPow2(combinedItems);

        while (combinedItems > newSize >> 1)
            newSize <<= 1;

        m_arraySize = newSize;

        auto oldKeys = std::move(m_keys);
        auto oldValues = std::move(m_values);
        m_items.store(0, std::memory_order_relaxed);

        try {
            m_keys.reset(new std::atomic<K>[m_arraySize]()); // zero init
            m_values.reset(new std::atomic<V>[m_arraySize]); // random init
        } catch (std::bad_alloc &e) {
            std::cerr << e.what() << std::endl;
            raise(SIGINT);
        }
        VTUNE_END_TASK(MultiMapMergeAllocate);


        VTUNE_TASK(MultiMapMergeCopy);
        tbb::parallel_for(tbb::blocked_range<std::size_t>(0, std::max(oldSize, other.capacity())),
                          [&oldKeys, &oldValues, oldSize, &other, &filter, this](const auto &r) {

                              K key = 0;
                              std::pair<K, V> pair;
                              for (auto i = r.begin(); i != r.end(); ++i) {
                                  if (i < oldSize) {
                                      key = oldKeys[i].load(std::memory_order_relaxed);

                                      if (key != 0 && !filter.count(key))
                                          this->_insert(key, oldValues[i].load(std::memory_order_relaxed));
                                  }

                                  if (i < other.capacity()) {
                                      pair = other.at(i);

                                      if (pair.first != 0 && !filter.count(pair.first))
                                          this->_insert(pair);
                                  }
                              }
                          });
    }

    const_iterator begin() const {
        return const_iterator(*this, 0, true);
    }

    const_iterator end() const {
        return const_iterator(*this, m_arraySize, false);
    }

    const_range_type range() const {
        return const_range_type(*this);
    }

private:

    void _insert(const K &key, const V &value) {
        ASSERT(key != 0);

        for (K idx = m_hasher(key); ; idx++) {

            idx &= m_arraySize - 1;
            ASSERT(idx < m_arraySize);

            // Load the key that was there.
            K probedKey = m_keys[idx].load();

            if (probedKey != 0) {
                continue; // Usually, it contains another key. Keep probing.
            }
            // The entry was free. Now let's try to take it using a CAS.
            K prevKey = 0;
            bool cas = m_keys[idx].compare_exchange_strong(prevKey, key);
            if (cas) {
                m_values[idx].store(value, std::memory_order_relaxed); // insert value
                ++m_items;
                return; // we just added the key to the set
            } else {
                continue; // another thread inserted in this position
            }
        }
    }

    void _insert(const std::pair<K, V> &pair) {
        insert(pair.first, pair.second);
    }

    void _migrate(const std::size_t idx, Concurrent_LP_MultiMap<K, V> &target) {
        auto pair = at(idx);
        if (pair.first != 0)
            target._insert(pair);
    }

private:
    std::size_t m_arraySize;
    std::atomic<std::size_t> m_items;
    std::unique_ptr<std::atomic<K>[]> m_keys;
    std::unique_ptr<std::atomic<V>[]> m_values;
    uint m_version;
    Hasher<K> m_hasher;

    std::atomic<std::size_t> m_currentCopyBlock;

    const static uint c_growThreshold = 10;

};

// conversion ctors

template<typename K, typename V>
LP_MultiMap<K, V>::LP_MultiMap(Concurrent_LP_MultiMap<K, V> &&other)
        : m_arraySize(other.m_arraySize),
          m_items(other.m_items),
          m_keys(reinterpret_cast<K *>(other.m_keys.release())),
          m_values(reinterpret_cast<V *>(other.m_values.release())),
          m_hasher(std::move(other.m_hasher)) { }

template<typename K, typename V>
Concurrent_LP_MultiMap<K, V>::Concurrent_LP_MultiMap(LP_MultiMap<K, V> &&other)
        : m_arraySize(other.m_arraySize),
          m_items(other.m_items),
          m_keys(reinterpret_cast<std::atomic<K> *>(other.m_keys.release())),
          m_values(reinterpret_cast<std::atomic<V> *>(other.m_values.release())),
          m_version(0),
          m_hasher(std::move(other.m_hasher)),
          m_currentCopyBlock(0) { }