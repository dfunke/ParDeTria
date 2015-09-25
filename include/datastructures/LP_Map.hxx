#pragma once

#include "LP_Set.hxx"

//forward declare concurrent version
template<typename K, typename V>
class Concurrent_LP_Map;

template<typename K, typename V>
class LP_Map {

    friend class Concurrent_LP_Map<K, V>;

public:

    LP_Map(std::size_t size, Hasher<K> hasher = Hasher<K>())
            : m_items(0),
              m_keys(nullptr),
              m_hasher(hasher) {
        // Initialize cells
        m_arraySize = nextPow2(size);

        m_keys.reset(new K[m_arraySize]()); // zero init
        m_values.reset(new V[m_arraySize]); // random init
    }

    LP_Map(LP_Map<K, V> &&other)
            : m_arraySize(other.m_arraySize),
              m_items(other.m_items),
              m_keys(std::move(other.m_keys)),
              m_values(std::move(other.m_values)),
              m_hasher(std::move(other.m_hasher)) { }

    //Conversion
    LP_Map(Concurrent_LP_Map<K, V> &&other);

    LP_Map<K, V> &operator=(LP_Map<K, V> &&other) {
        m_arraySize = other.m_arraySize;
        m_items = other.m_items;
        m_keys = std::move(other.m_keys);
        m_values = std::move(other.m_values);
        m_hasher = std::move(other.m_hasher);

        return *this;
    }

    bool insert(const K &key, const V &value) {
        ASSERT(key != 0);

        if (!m_rehashing && m_items / m_arraySize > 0.5)
            rehash(m_arraySize << 1);

        for (K idx = m_hasher(key); ; idx++) {
            idx &= m_arraySize - 1;
            ASSERT(idx < m_arraySize);

            // Load the key that was there.
            K probedKey = m_keys[idx];

            if (probedKey == key) {
                m_values[idx] = value; // key already exists, update value
                return false; // we didn't insert a new key
            }
            else {
                // The entry was either free, or contains another key.
                if (probedKey != 0)
                    continue; // Usually, it contains another key. Keep probing.
                // The entry was free. take it
                m_keys[idx] = key;
                m_values[idx] = value;
                ++m_items;
                return true;
            }
        }
    }

    bool insert(const std::pair<K, V> &pair) {
        return insert(pair.first, pair.second);
    }

    V get(const K &key) const {
        ASSERT(key != 0);

        for (K idx = m_hasher(key); ; idx++) {
            idx &= m_arraySize - 1;
            K probedKey = m_keys[idx];
            if (probedKey == key)
                return m_values[idx];;
            if (probedKey == 0)
                return 0;
        }

    }

    bool contains(const K &key) const {
        return get(key) != 0;
    }

    std::size_t count(const K &key) const {
        return contains(key);
    }


    bool empty() const { return m_items == 0; };

    bool empty(const std::size_t idx) const {
        return m_keys[idx] == 0;
    };

    std::pair<K, V> at(const std::size_t idx) const {

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

    void merge(LP_Map<K, V> &&other) {

        rehash((capacity() + other.capacity()) << 1);

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
    void merge(LP_Map<K, V> &&other, const Set &filter) {
        rehash((capacity() + other.capacity()) << 1, filter);

        for (std::size_t i = 0; i < other.capacity(); ++i) {
            if (other.m_keys[i] != 0 && !filter.count(other.m_keys[i]))
                insert(other.m_keys[i], other.m_values[i]);
        }
    }


public:
    typedef _detail::iterator<const LP_Map<K, V>, std::pair<const K, V>> const_iterator;
    typedef _detail::range_type<const LP_Map<K, V>, const_iterator> const_range_type;

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
class Concurrent_LP_Map {

    friend class LP_Map<K, V>;

    friend class GrowingHashTable<Concurrent_LP_Map<K, V>>;

    friend class GrowingHashTableHandle<Concurrent_LP_Map<K, V>>;

    friend class ConstGrowingHashTableHandle<Concurrent_LP_Map<K, V>>;

public:

    Concurrent_LP_Map(std::size_t size, const uint version = 0,
                      Hasher<K> hasher = Hasher<K>())
            : m_items(0),
              m_keys(nullptr),
              m_version(version),
              m_hasher(hasher),
              m_currentCopyBlock(0) {
        // Initialize cells
        m_arraySize = nextPow2(size);
        m_keys.reset(new std::atomic<K>[m_arraySize]()); // zero init
        m_values.reset(new std::atomic<V>[m_arraySize]); // random init
    }

    Concurrent_LP_Map(Concurrent_LP_Map<K, V> &&other)
            : m_arraySize(other.m_arraySize),
              m_items(other.m_items.load()),
              m_keys(std::move(other.m_keys)),
              m_values(std::move(other.m_values)),
              m_version(other.m_version),
              m_hasher(std::move(other.m_hasher)),
              m_currentCopyBlock(0) { }

    //conversion
    Concurrent_LP_Map(LP_Map<K, V> &&other);

    Concurrent_LP_Map<K, V> &operator=(Concurrent_LP_Map<K, V> &&other) {
        m_arraySize = other.m_arraySize;
        m_items.store(other.m_items.load());
        m_keys = std::move(other.m_keys);
        m_values = std::move(other.m_values);
        m_version = other.m_version;
        m_hasher = std::move(other.m_hasher);
        m_currentCopyBlock.store(other.m_currentCopyBlock.load());

        return *this;
    }


    InsertReturn insert(const K &key, const V &value) {
        ASSERT(key != 0);

        if (m_items.load() > m_arraySize >> 1)
            return InsertReturn::State::Full;

        for (K idx = m_hasher(key); ; idx++) {
            idx &= m_arraySize - 1;
            ASSERT(idx < m_arraySize);

            // Load the key that was there.
            K probedKey = m_keys[idx].load();

            if (probedKey == key) {
                m_values[idx].store(value, std::memory_order_relaxed); // update value
                return false; // the key is already in the set, return false;
            }
            else {
                // The entry was either free, or contains another key.
                if (probedKey != 0)
                    continue; // Usually, it contains another key. Keep probing.
                // The entry was free. Now let's try to take it using a CAS.
                K prevKey = 0;
                bool cas = m_keys[idx].compare_exchange_strong(prevKey, key);
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

    InsertReturn insert(const std::pair<K, V> &pair) {
        return insert(pair.first, pair.second);
    }

    V get(const K &key) const {
        ASSERT(key != 0);
        for (K idx = m_hasher(key); ; idx++) {
            idx &= m_arraySize - 1;
            K probedKey = m_keys[idx].load();
            if (probedKey == key)
                return m_values[idx].load();;
            if (probedKey == 0)
                return 0;
        }

    }

    bool contains(const K &key) const {
        return get(key) != 0;
    }

    std::size_t count(const K &key) const {
        return contains(key);
    }

    bool empty() const { return m_items.load() == 0; };

    bool empty(const std::size_t idx) const { return m_keys[idx].load(std::memory_order_relaxed) == 0; };

    std::pair<K, V> at(const std::size_t idx) const {
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

        m_keys.reset(new std::atomic<K>[m_arraySize]()); //zero init
        m_values.reset(new std::atomic<V>[m_arraySize]); //random init

        tbb::parallel_for(std::size_t(0), oldSize, [&oldKeys, &oldValues, this](const K i) {
            if (oldKeys[i].load(std::memory_order_relaxed) != 0)
                insert(oldKeys[i].load(std::memory_order_relaxed), oldValues[i].load(std::memory_order_relaxed));
        });
    }

    void unsafe_merge(Concurrent_LP_Map<K, V> &&other) {

        unsafe_rehash((capacity() + other.capacity()) << 1);

        tbb::parallel_for(std::size_t(0), other.capacity(), [&other, this](const K i) {
            if (other.m_keys[i].load(std::memory_order_relaxed) != 0)
                insert(other.m_keys[i].load(std::memory_order_relaxed),
                       other.m_values[i].load(std::memory_order_relaxed));
        });
    }


    template<class Set>
    void unsafe_merge(Concurrent_LP_Map<K, V> &&other, const Set &filter) {

        std::size_t oldSize = m_arraySize;
        m_arraySize = nextPow2((capacity() + other.capacity()) << 1);

        auto oldKeys = std::move(m_keys);
        auto oldValues = std::move(m_values);
        m_items.store(0, std::memory_order_relaxed);

        m_keys.reset(new std::atomic<K>[m_arraySize]()); // zero init
        m_values.reset(new std::atomic<V>[m_arraySize]); // random init

        tbb::parallel_for(tbb::blocked_range<std::size_t>(0, std::max(oldSize, other.capacity())),
                          [&oldKeys, &oldValues, oldSize, &other, &filter, this](const auto &r) {

                              K val = 0;
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
    typedef _detail::iterator<const Concurrent_LP_Map<K, V>, std::pair<const K, V>> const_iterator;
    typedef _detail::range_type<const Concurrent_LP_Map<K, V>, const_iterator> const_range_type;

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

            if (probedKey == key) {
                m_values[idx].store(value, std::memory_order_relaxed); // update value
            }
            else {
                // The entry was either free, or contains another key.
                if (probedKey != 0)
                    continue; // Usually, it contains another key. Keep probing.
                // The entry was free. Now let's try to take it using a CAS.
                K prevKey = 0;
                bool cas = m_keys[idx].compare_exchange_strong(prevKey, key);
                if (cas) {
                    m_values[idx].store(value, std::memory_order_relaxed); // insert value
                    ++m_items;
                }
                else if (prevKey == key) {
                    m_values[idx].store(value, std::memory_order_relaxed); // update value
                } else
                    continue; // another thread inserted a different key in this position
            }
        }
    }

    void _insert(const std::pair<K, V> &pair) {
        insert(pair.first, pair.second);
    }

    void _migrate(const std::size_t idx, Concurrent_LP_Map<K, V> &target) {
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

};

//conversion constructors

template<typename K, typename V>
LP_Map<K, V>::LP_Map(Concurrent_LP_Map<K, V> &&other)
        : m_arraySize(other.m_arraySize),
          m_items(other.m_items),
          m_keys(reinterpret_cast<K *>(other.m_keys.release())),
          m_values(reinterpret_cast<V *>(other.m_values.release())),
          m_hasher(std::move(other.m_hasher)) { }

template<typename K, typename V>
Concurrent_LP_Map<K, V>::Concurrent_LP_Map(LP_Map<K, V> &&other)
        : m_arraySize(other.m_arraySize),
          m_items(other.m_items),
          m_keys(reinterpret_cast<std::atomic<K> *>(other.m_keys.release())),
          m_values(reinterpret_cast<std::atomic<V> *>(other.m_values.release())),
          m_version(0),
          m_hasher(std::move(other.m_hasher)),
          m_currentCopyBlock(0) { }

