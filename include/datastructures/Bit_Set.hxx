#pragma once

#include "datastructures/LP_Set.hxx"

#include "utils/Misc.h"
#include "utils/ASSERT.h"

namespace _detail {

    template<class Container, typename Value>
    struct bitset_iterator : public std::iterator<std::bidirectional_iterator_tag, Value> {

    private:
        static const uint cENTRIES = (CHAR_BIT * sizeof(typename Container::tKeyType));

    public:
        bitset_iterator(Container &container, std::size_t idx, bool findNext)
                : m_container(container),
                  m_idx(m_container._idx(idx)),
                  m_lowerBound(container.lowerBound()),
                  m_end(m_container.upperBound() - m_container.lowerBound()) {

            if (findNext)
                _advance(true);

        }

        bitset_iterator(const bitset_iterator &o)
                : m_container(o.m_container),
                  m_idx(o.m_idx),
                  m_lowerBound(o.m_lowerBound),
                  m_end(o.m_end) { }

        bitset_iterator &operator++() {
            _advance(false);
            return *this;
        }

        bitset_iterator operator++(int) {
            auto tmp = *this;
            ++(*this);
            return tmp;
        }

        bitset_iterator &operator--() {
            _dec(false);
            return *this;
        }

        bitset_iterator operator--(int) {
            auto tmp = *this;
            --(*this);
            return tmp;
        }

        std::size_t operator-(const bitset_iterator &other) const {
            return m_idx - other.m_idx;
        }

        std::size_t operator+(const bitset_iterator &other) const {
            return m_idx + other.m_idx;
        }

        std::size_t operator+(const std::size_t &a) const {
            return m_idx + a;
        }

        bitset_iterator operator=(const bitset_iterator &other) {
            m_idx = other.m_idx;
            return *this;
        }

        bool operator==(const bitset_iterator &j) const {
            return m_idx == j.m_idx;
        }

        bool operator!=(const bitset_iterator &j) const { return !(*this == j); }

        bool operator<(const bitset_iterator &other) const {
            return m_idx < other.m_idx;
        }

        bool operator<=(const bitset_iterator &other) const {
            return m_idx <= other.m_idx;
        }

        bool operator>(const bitset_iterator &other) const {
            return m_idx > other.m_idx;
        }

        bool operator>=(const bitset_iterator &other) const {
            return m_idx >= other.m_idx;
        }

        auto operator*() const { return m_lowerBound + m_idx; }

        //const auto *operator->() const { return &m_container.at(m_idx); }

        void half(bitset_iterator &begin, bitset_iterator & end){
            _setIdx(begin + ((end - begin) / 2), true);
        }

    private:

        void _setIdx(const std::size_t &idx, bool findNext) {
            m_idx = idx;
            if (findNext)
                _advance(true);
        }

        void _advance(const bool testFirst) {
            if (testFirst) {
                while ((m_idx < m_end && !m_container._isSet(m_idx)))
                    m_idx += m_idx % cENTRIES != 0 ?
                             1 // we are in the middle of a byte
                                                   : m_container.m_array[m_idx / cENTRIES] ? 1
                                                                                           : cENTRIES; // maybe we can skip a byte
            } else {
                while (m_idx += m_idx % cENTRIES != 0 ?
                                1 // we are in the middle of a byte
                                                      : m_container.m_array[m_idx / cENTRIES] ? 1
                                                                                              : cENTRIES, // maybe we can skip a byte
                        (m_idx < m_end && !m_container._isSet(m_idx))) { }
            }
        }

        void _dec(const bool testFirst) {
            if (testFirst) {
                while ((m_idx >= 0 && !m_container._isSet(m_idx)))
                    --m_idx;
            } else {
                while (--m_idx, (m_idx >= 0 && !m_container._isSet(m_idx))) { }
            }
        }

    protected:
        Container &m_container;
        std::size_t m_idx;
        const std::size_t m_lowerBound;
        const std::size_t m_end;
    };
}

class Concurrent_Bit_Set;

class Bit_Set {

    friend class Concurrent_Bit_Set;

public:
    typedef uint tKeyType;

public:

    Bit_Set(std::size_t lowerBound, std::size_t upperBound, bool zero = true)
            : m_array(nullptr),
              m_ones(0) {

        // Initialize cells
        m_lowerBound = (lowerBound / cENTRIES) * cENTRIES; // round down
        m_upperBound = ((upperBound + cENTRIES - 1) / cENTRIES) * cENTRIES; // round up

        m_arraySize = (m_upperBound - m_lowerBound) / cENTRIES;
        try {
            if (zero)
                m_array.reset(new tKeyType[m_arraySize]()); //value init
            else
                m_array.reset(new tKeyType[m_arraySize]); //random init
        } catch (std::bad_alloc &e) {
            std::cerr << e.what() << std::endl;
            raise(SIGINT);
        }
    }

    Bit_Set(Bit_Set &&other)
            : m_arraySize(other.m_arraySize),
              m_lowerBound(other.m_lowerBound),
              m_upperBound(other.m_upperBound),
              m_array(std::move(other.m_array)),
              m_ones(other.m_ones) { }

    //Conversion
    Bit_Set(Concurrent_Bit_Set &&other);

    Bit_Set &operator=(Bit_Set &&other) {
        m_arraySize = other.m_arraySize;
        m_lowerBound = other.m_lowerBound;
        m_upperBound = other.m_upperBound;
        m_array = std::move(other.m_array);
        m_ones = other.m_ones;

        return *this;
    }

    bool testAndSet(const tKeyType &idx) {
        RAISE(m_lowerBound <= idx && idx < m_upperBound);
        RAISE(_block(idx) < m_arraySize);

        bool test = isSet(idx);
        m_array[_block(idx)] |= 1 << (idx % cENTRIES);
        m_ones += !test;

        return test;
    }

    bool testAndUnset(const tKeyType &idx) {
        RAISE(m_lowerBound <= idx && idx < m_upperBound);
        RAISE(_block(idx) < m_arraySize);

        bool test = isSet(idx);
        m_array[_block(idx)] &= ~(1 << (idx % cENTRIES));
        m_ones -= test;

        return test;
    }

    void batchSet(const tKeyType &from, const tKeyType &to) {
        RAISE(m_lowerBound <= from && from < m_upperBound);
        RAISE(m_lowerBound <= to && to <= m_upperBound);

        /*auto print = [this] (const std::string & msg) {
            std::cout << msg << ": ";
            for(std::size_t i = 0; i < capacity(); ++i)
                std::cout << isSet(i);
            std::cout << std::endl;
        };*/

        std::size_t fromBlock = _block(from);
        std::size_t toBlock = _block(to - 1); // to not included

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

        if (m_lowerBound <= idx && idx < m_upperBound) {
            RAISE(_block(idx) < m_arraySize);

            return m_array[_block(idx)] && (m_array[_block(idx)] & (1 << (idx % cENTRIES)));
        } else {
            return false;
        }
    }

    bool zero() const { return m_ones == 0; }

    bool one() const { return m_ones == m_upperBound - m_lowerBound; }

    bool zero(const std::size_t idx) const {
        return !isSet(idx);
    }

    bool one(const std::size_t idx) const {
        return isSet(idx);
    }

    //std::size_t capacity() const { return m_upperBound - m_lowerBound; }
    std::size_t lowerBound() const { return m_lowerBound; }

    std::size_t upperBound() const { return m_upperBound; }

    std::size_t zeros() const { return m_upperBound - m_lowerBound - m_ones; }

    std::size_t ones() const { return m_ones; }

    template<class BitSet>
    void copy(const BitSet &other) {

        m_arraySize = other.m_arraySize;
        m_lowerBound = other.m_lowerBound;
        m_upperBound = other.m_upperBound;
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

    void resize(const std::size_t &upperBound) {
        resize(m_lowerBound, upperBound);
    }

    void resize(const std::size_t &lowerBound, const std::size_t &upperBound) {

        tArrayPtr oldArray = std::move(m_array);
        std::size_t oldLowerBound = m_lowerBound;
        std::size_t oldUpperBound = m_upperBound;

        m_lowerBound = std::min(m_lowerBound, (lowerBound / cENTRIES) * cENTRIES); // round down
        m_upperBound = std::max(m_upperBound, ((upperBound + cENTRIES - 1) / cENTRIES) * cENTRIES); // round up

        oldLowerBound = _block(oldLowerBound);
        oldUpperBound = _block(oldUpperBound);

        m_arraySize = (m_upperBound - m_lowerBound) / cENTRIES;

        try {
            m_array.reset(new tKeyType[m_arraySize]); //random init
        } catch (std::bad_alloc &e) {
            std::cerr << e.what() << std::endl;
            raise(SIGINT);
        }

        std::size_t i;
        for (i = 0; i < oldLowerBound; ++i) {
            m_array[i] = (tKeyType) 0;
        }
        for (; i < oldUpperBound; ++i) {
            m_array[i] = oldArray[i - oldLowerBound];
        }
        for (; i < m_arraySize; ++i) {
            m_array[i] = (tKeyType) 0;
        }
    }

    template<class BitSet>
    void merge(BitSet &&other) {

        tArrayPtr oldArray = std::move(m_array);
        std::size_t oldLowerBound = m_lowerBound;
        std::size_t oldUpperBound = m_upperBound;

        m_lowerBound = std::min(oldLowerBound, other.m_lowerBound);
        m_upperBound = std::max(oldUpperBound, other.m_upperBound);

        m_arraySize = (m_upperBound - m_lowerBound) / cENTRIES;

        try {
            m_array.reset(new tKeyType[m_arraySize]); //random init
        } catch (std::bad_alloc &e) {
            std::cerr << e.what() << std::endl;
            raise(SIGINT);
        }

        m_ones = 0;
        for (std::size_t i = m_lowerBound;
             i < m_upperBound;
             i += cENTRIES) {
            m_array[_block(i)] = ((oldLowerBound <= i && i < oldUpperBound) ? oldArray[(i - oldLowerBound) / cENTRIES]
                                                                            : (tKeyType) 0)
                                 | ((other.m_lowerBound <= i && i < other.m_upperBound) ? other.m_array[other._block(i)]
                                                                                        : (tKeyType) 0);
            m_ones += popcount(m_array[_block(i)]);
        }

    }

    template<class BitSet>
    void filter(const BitSet &filter) {

        m_ones = 0;
        for (std::size_t i = std::max(m_lowerBound, filter.m_lowerBound);
             i < std::min(m_upperBound, filter.m_upperBound);
             i += cENTRIES) {
            m_array[_block(i)] &= ~filter.m_array[filter._block(i)];
            m_ones += popcount(m_array[_block(i)]);
        }
    }

    template<class BitSet, class Filter>
    void mergeFilter(BitSet &&other, const Filter &filter) {

        tArrayPtr oldArray = std::move(m_array);
        std::size_t oldLowerBound = m_lowerBound;
        std::size_t oldUpperBound = m_upperBound;

        m_lowerBound = std::min(oldLowerBound, other.m_lowerBound);
        m_upperBound = std::max(oldUpperBound, other.m_upperBound);

        m_arraySize = (m_upperBound - m_lowerBound) / cENTRIES;

        try {
            m_array.reset(new tKeyType[m_arraySize]); //random init
        } catch (std::bad_alloc &e) {
            std::cerr << e.what() << std::endl;
            raise(SIGINT);
        }

        m_ones = 0;
        for (std::size_t i = m_lowerBound;
             i < m_upperBound;
             i += cENTRIES) {
            m_array[_block(i)] = (((oldLowerBound <= i && i < oldUpperBound) ? oldArray[(i - oldLowerBound) / cENTRIES]
                                                                             : (tKeyType) 0)
                                  |
                                  ((other.m_lowerBound <= i && i < other.m_upperBound) ? other.m_array[other._block(i)]
                                                                                       : (tKeyType) 0))
                                 &
                                 ~(filter.m_lowerBound <= i && i < filter.m_upperBound ? filter.m_array[filter._block(
                                         i)] : (tKeyType) 0);
            m_ones += popcount(m_array[_block(i)]);
        }
    }



    // compat definitions
    bool insert(const tKeyType &key) { return !testAndSet(key); }

    bool erase(const tKeyType &key) { return testAndUnset(key); }

    bool contains(const tKeyType &key) const { return isSet(key); }

    std::size_t count(const tKeyType &key) const { return contains(key); }

    bool empty() const { return zero(); };

    bool empty(const std::size_t idx) const { return zero(idx); }

    tKeyType at(const std::size_t idx) const { return idx * isSet(idx); }

    std::size_t size() const { return ones(); }

    std::size_t capacity() const { return upperBound(); }

private:
    inline std::size_t _idx(const std::size_t &idx) const {
        return idx - m_lowerBound;
    }

    inline std::size_t _block(const std::size_t &idx) const {
        return _idx(idx) / cENTRIES;
    }

    bool _isSet(const tKeyType &idx) const {
        return m_array[idx / cENTRIES] && (m_array[idx / cENTRIES] & (1 << (idx % cENTRIES)));
    }

public:
    typedef _detail::bitset_iterator<const Bit_Set, tKeyType> const_iterator;

    friend const_iterator;

    typedef _detail::range_type<const Bit_Set, const_iterator> const_range_type;

    const_iterator begin() const {
        return const_iterator(*this, m_lowerBound, true);
    }

    const_iterator end() const {
        return const_iterator(*this, m_upperBound, false);
    }

    const_range_type range() const {
        return const_range_type(*this);
    }

private:
    typedef std::unique_ptr<tKeyType[]> tArrayPtr;

private:
    std::size_t m_arraySize;

    std::size_t m_lowerBound;
    std::size_t m_upperBound;

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

    Concurrent_Bit_Set(std::size_t lowerBound, std::size_t upperBound, bool zero = true)
            : m_array(nullptr),
              m_ones(0) {

        // Initialize cells
        m_lowerBound = (lowerBound / cENTRIES) * cENTRIES; // round down
        m_upperBound = ((upperBound + cENTRIES - 1) / cENTRIES) * cENTRIES; // round up

        m_arraySize = (m_upperBound - m_lowerBound) / cENTRIES;
        try {
            if (zero)
                m_array.reset(new std::atomic<tKeyType>[m_arraySize]()); //value init
            else
                m_array.reset(new std::atomic<tKeyType>[m_arraySize]); //random init
        } catch (std::bad_alloc &e) {
            std::cerr << e.what() << std::endl;
            raise(SIGINT);
        }
    }

    Concurrent_Bit_Set(Concurrent_Bit_Set &&other)
            : m_arraySize(other.m_arraySize),
              m_lowerBound(other.m_lowerBound),
              m_upperBound(other.m_upperBound),
              m_array(std::move(other.m_array)),
              m_ones(other.m_ones.load()) { }

    //Conversion
    Concurrent_Bit_Set(Bit_Set &&other);

    Concurrent_Bit_Set &operator=(Concurrent_Bit_Set &&other) {
        m_arraySize = other.m_arraySize;
        m_lowerBound = other.m_lowerBound;
        m_upperBound = other.m_upperBound;
        m_array = std::move(other.m_array);
        m_ones.store(other.m_ones.load());

        return *this;
    }

    bool testAndSet(const tKeyType &idx) {
        RAISE(m_lowerBound <= idx && idx < m_upperBound);
        RAISE(_block(idx) < m_arraySize);

        std::size_t i = _block(idx);
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

    bool testAndUnset(const tKeyType &idx) {
        RAISE(m_lowerBound <= idx && idx < m_upperBound);
        RAISE(_block(idx) < m_arraySize);

        std::size_t i = _block(idx);
        tKeyType mask = (1 << (idx % cENTRIES));

        tKeyType des, val = m_array[i].load();

        while ((val & mask)) {
            des = val & ~mask;

            bool cas = m_array[i].compare_exchange_strong(val, des);

            if (cas) {
                --m_ones; // we just unsset the bit, dec m_ones
                return true; // the bit was set before
            }
        }

        // value was already zero
        return false;
    }

    bool isSet(const tKeyType &idx) const {

        if (m_lowerBound <= idx && idx < m_upperBound) {
            RAISE(_block(idx) < m_arraySize);

            tKeyType val = m_array[_block(idx)].load();
            return val && (val & (1 << (idx % cENTRIES)));
        } else {
            return false;
        }
    }

    bool zero() const { return m_ones.load() == 0; }

    bool one() const { return m_ones.load() == m_upperBound - m_lowerBound; }

    bool zero(const std::size_t idx) const {
        return !isSet(idx);
    }

    bool one(const std::size_t idx) const {
        return isSet(idx);
    }

    //std::size_t capacity() const { return m_upperBound - m_lowerBound; }
    std::size_t lowerBound() const { return m_lowerBound; }

    std::size_t upperBound() const { return m_upperBound; }

    std::size_t zeros() const { return m_upperBound - m_lowerBound - m_ones.load(); }

    std::size_t ones() const { return m_ones.load(); }

    // compat definitions
    bool insert(const tKeyType &key) { return !testAndSet(key); }

    bool erase(const tKeyType &key) { return testAndUnset(key); }

    bool contains(const tKeyType &key) const { return isSet(key); }

    std::size_t count(const tKeyType &key) const { return contains(key); }

    bool empty() const { return zero(); };

    bool empty(const std::size_t idx) const { return zero(idx); }

    tKeyType at(const std::size_t idx) const { return idx * isSet(idx); }

    std::size_t size() const { return ones(); }

    std::size_t capacity() const { return upperBound(); }

private:
    inline std::size_t _idx(const std::size_t &idx) const {
        return idx - m_lowerBound;
    }

    inline std::size_t _block(const std::size_t &idx) const {
        return _idx(idx) / cENTRIES;
    }

    bool _isSet(const tKeyType &idx) const {
        tKeyType val = m_array[idx / cENTRIES].load();
        return val && (val & (1 << (idx % cENTRIES)));
    }

public:
    typedef _detail::bitset_iterator<const Concurrent_Bit_Set, tKeyType> const_iterator;

    friend const_iterator;

    typedef _detail::range_type<const Concurrent_Bit_Set, const_iterator> const_range_type;

    const_iterator begin() const {
        return const_iterator(*this, m_lowerBound, true);
    }

    const_iterator end() const {
        return const_iterator(*this, m_upperBound, false);
    }

    const_range_type range() const {
        return const_range_type(*this);
    }

private:
    typedef std::unique_ptr<std::atomic<tKeyType>[]> tArrayPtr;

private:
    std::size_t m_arraySize;
    std::size_t m_lowerBound;
    std::size_t m_upperBound;

    tArrayPtr m_array;

    std::atomic<std::size_t> m_ones;

private:
    static const uint cENTRIES = (CHAR_BIT * sizeof(tKeyType));
};

// conversion ctors


Bit_Set::Bit_Set(Concurrent_Bit_Set &&other)
        : m_arraySize(other.m_arraySize),
          m_lowerBound(other.m_lowerBound),
          m_upperBound(other.m_upperBound),
          m_array(reinterpret_cast<uint *>(other.m_array.release())),
          m_ones(other.m_ones.load()) { }

Concurrent_Bit_Set::Concurrent_Bit_Set(Bit_Set &&other)
        : m_arraySize(other.m_arraySize),
          m_lowerBound(other.m_lowerBound),
          m_upperBound(other.m_upperBound),
          m_array(reinterpret_cast<std::atomic<uint> *>(other.m_array.release())),
          m_ones(other.m_ones) { }
