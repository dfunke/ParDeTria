#pragma once

#include <memory>
#include <vector>
#include <algorithm>

#include "datastructures/LP_Set.hxx"

template<typename T, typename IDX>
class BlockedArray2;

namespace _detail {

    template<class Container, typename Value>
    struct block_iterator : public std::iterator<std::bidirectional_iterator_tag, Value> {

    public:
        block_iterator(Container &container, std::size_t block, std::size_t idx)
                : m_container(container),
                  m_idx(idx),
                  m_block(block) {
        }

        block_iterator(const block_iterator &o)
                : m_container(o.m_container),
                  m_idx(o.m_idx),
                  m_block(o.m_block) { }

        block_iterator &operator++() {
            _advance();
            return *this;
        }

        block_iterator operator++(int) {
            auto tmp = *this;
            ++(*this);
            return tmp;
        }

        block_iterator &operator--() {
            _dec();
            return *this;
        }

        block_iterator operator--(int) {
            auto tmp = *this;
            --(*this);
            return tmp;
        }

        std::size_t operator-(const block_iterator &other) const {
            return m_idx - other.m_idx;
        }

        std::size_t operator+(const block_iterator &other) const {
            return m_idx + other.m_idx;
        }

        block_iterator operator=(const block_iterator &other) {
            m_idx = other.m_idx;
            m_block = other.m_block;
            return *this;
        }

        bool operator==(const block_iterator &j) const {
            return m_block == j.m_block && m_idx == j.m_idx;
        }

        bool operator!=(const block_iterator &j) const { return !(*this == j); }

        bool operator<(const block_iterator &j) const {
            return m_block < j.m_block || (m_block == j.m_block && m_idx < j.m_idx);
        }

        bool operator<=(const block_iterator &j) const {
            return m_block < j.m_block || (m_block == j.m_block && m_idx <= j.m_idx);
        }

        bool operator>(const block_iterator &j) const {
            return m_block > j.m_block || (m_block == j.m_block && m_idx > j.m_idx);;
        }

        bool operator>=(const block_iterator &j) const {
            return m_block > j.m_block || (m_block == j.m_block && m_idx >= j.m_idx);;
        }

        auto &operator*() const { return m_container.m_blocks[m_block]->at(m_idx); }

        auto *operator->() const { return &m_container.m_blocks[m_block]->at(m_idx); }

        void half(block_iterator &begin, block_iterator & end){
            _setIdx((end + begin) / 2, true);
        }

    private:

        void _setIdx(const std::size_t &idx, __attribute__((unused)) bool findNext) {
            m_block = m_container._findBlock(idx);
            m_idx = std::max(idx, m_container.m_blocks[m_block]->min());
        }

        void _advance() {
            ++m_idx;
            if (m_container.m_blocks[m_block]->max() == m_idx && m_block < m_container.m_blocks.size() - 1) {
                ++m_block;
                m_idx = m_container.m_blocks[m_block]->min();
            }
        }

        void _dec() {
            if (m_container.m_blocks[m_block]->min() == m_idx && m_block > 0) {
                --m_block;
                m_idx = m_container.m_blocks[m_block]->max();
            }
            --m_idx;
        }

    protected:
        Container &m_container;
        std::size_t m_idx;
        std::size_t m_block;
    };

    template<class Container, typename Value>
    struct filtered_block_iterator : public std::iterator<std::bidirectional_iterator_tag, Value> {

    public:
        filtered_block_iterator(Container &container, std::size_t block, std::size_t idx, const bool findNext)
                : m_container(container),
                  m_idx(idx),
                  m_block(block) {

            if (findNext)
                _advance(true);
        }

        filtered_block_iterator(const filtered_block_iterator &o)
                : m_container(o.m_container),
                  m_idx(o.m_idx),
                  m_block(o.m_block) { }

        filtered_block_iterator &operator++() {
            _advance(false);
            return *this;
        }

        filtered_block_iterator operator++(int) {
            auto tmp = *this;
            ++(*this);
            return tmp;
        }

        filtered_block_iterator &operator--() {
            _dec(false);
            return *this;
        }

        filtered_block_iterator operator--(int) {
            auto tmp = *this;
            --(*this);
            return tmp;
        }

        std::size_t operator-(const filtered_block_iterator &other) const {
            return m_idx - other.m_idx;
        }

        std::size_t operator+(const filtered_block_iterator &other) const {
            return m_idx + other.m_idx;
        }

        filtered_block_iterator operator=(const filtered_block_iterator &other) {
            m_idx = other.m_idx;
            m_block = other.m_block;
            return *this;
        }

        bool operator==(const filtered_block_iterator &j) const {
            return m_block == j.m_block && m_idx == j.m_idx;
        }

        bool operator!=(const filtered_block_iterator &j) const { return !(*this == j); }

        bool operator<(const filtered_block_iterator &j) const {
            return m_block < j.m_block || (m_block == j.m_block && m_idx < j.m_idx);
        }

        bool operator<=(const filtered_block_iterator &j) const {
            return m_block < j.m_block || (m_block == j.m_block && m_idx <= j.m_idx);
        }

        bool operator>(const filtered_block_iterator &j) const {
            return m_block > j.m_block || (m_block == j.m_block && m_idx > j.m_idx);;
        }

        bool operator>=(const filtered_block_iterator &j) const {
            return m_block > j.m_block || (m_block == j.m_block && m_idx >= j.m_idx);;
        }

        auto &operator*() const { return m_container.m_blocks[m_block]->at(m_idx); }

        auto *operator->() const { return &m_container.m_blocks[m_block]->at(m_idx); }

        void half(filtered_block_iterator &begin, filtered_block_iterator & end){
            _setIdx((end + begin) / 2, true);
        }

    private:

        void _setIdx(const std::size_t &idx, bool findNext) {
            m_block = m_container._findBlock(idx);
            m_idx = std::max(idx, std::size_t(m_container.m_blocks[m_block]->min()));
            if (findNext)
                _advance(true);
        }

        inline void _advance(const bool testFirst) {
            if (!testFirst) { //advance if we don't want to test first
                _incIdx();
            }

            while (m_idx < m_container.m_blocks.back()->max()
                   && !m_container._valid(m_container.m_blocks[m_block]->at(m_idx))) {

                _incIdx();
            }
        }

        inline void _incIdx() {
            ++m_idx;
            if (m_container.m_blocks[m_block]->max() == m_idx && m_block < m_container.m_blocks.size() - 1) {
                ++m_block;
                m_idx = m_container.m_blocks[m_block]->min();
            }
        }

        inline void _dec(const bool testFirst) {
            if (!testFirst) { //decrement if we don't test first
                _decIdx();
            }

            while (m_idx > m_container.m_blocks.front()->min()
                   && !m_container._valid(m_container.m_blocks[m_block]->at(m_idx))) {

                _decIdx();
            }
        }

        inline void _decIdx() {
            if (m_container.m_blocks[m_block]->min() == m_idx && m_block > 0) {
                --m_block;
                m_idx = m_container.m_blocks[m_block]->max();
            }
            --m_idx;
        }

    protected:
        Container &m_container;
        std::size_t m_idx;
        std::size_t m_block;
    };

    template<typename T, typename IDX = std::size_t>
    class Block {

    private:
        typedef std::unique_ptr<T[]> data_ptr;

    public:
        Block(const IDX min, const IDX max)
                : m_min(min), m_max(max),
                  m_data(new T[max - min]()) { }

        Block(Block &&other) {
            m_min = other.m_min;
            m_max = other.m_max;
            m_data = std::move(other.m_data);
        }

        const T &operator[](const IDX idx) const {
            return at(idx);
        }

        T &operator[](const IDX idx) {
            return at(idx);
        }

        const T &at(const IDX idx) const {
#ifndef NDEBUG
            if (!(m_min <= idx && idx < m_max)) {
                throw std::out_of_range(
                        "index: " + std::to_string(idx) + " size: " + std::to_string(m_min) + " to " +
                        std::to_string(m_max));
            }
#endif
            return m_data[idx - m_min];
        }

        T &at(const IDX idx) {
#ifndef NDEBUG
            if (!(m_min <= idx && idx < m_max)) {
                throw std::out_of_range(
                        "index: " + std::to_string(idx) + " size: " + std::to_string(m_min) + " to " +
                        std::to_string(m_max));
            }
#endif
            return m_data[idx - m_min];
        }

        bool contains(const IDX idx) const {
            return m_min <= idx && idx < m_max;
        }

        IDX min() const {
            return m_min;
        }

        IDX max() const {
            return m_max;
        }

    private:
        friend class boost::serialization::access;

        template<class Archive>
        void save(Archive &ar, __attribute__((unused)) const unsigned int version) const {
            // invoke serialization of the base class
            ar << m_min;
            ar << m_max;

            ar << boost::serialization::make_array<T>(m_data.get(), size());
        }

        template<class Archive>
        void load(Archive &ar, __attribute__((unused)) const unsigned int version) {
            // invoke serialization of the base class
            ar >> m_min;
            ar >> m_max;

            ar >> boost::serialization::make_array<T>(m_data.get(), size());
        }

        BOOST_SERIALIZATION_SPLIT_MEMBER()


    private:
        IDX m_min;
        IDX m_max;

        data_ptr m_data;
    };
}


template<typename T, typename IDX = std::size_t>
class BlockedArray2 {

public:
    typedef _detail::Block<T, IDX> block;
    typedef std::unique_ptr<block> block_ptr;

public:
    BlockedArray2(const IDX min, const IDX max) {
        m_blocks.emplace_back(std::make_unique<block>(min, max));
    }

    BlockedArray2(BlockedArray2 &&other) {
        m_blocks = std::move(other.m_blocks);
    }

    BlockedArray2 &operator=(BlockedArray2 &&other) {
        m_blocks = std::move(other.m_blocks);
        m_hint = other.m_hint;

        return *this;
    }

    block &addBlock(const IDX min, const IDX max) {
        block_ptr nBlock = std::make_unique<block>(min, max);
        block &nBlockRef = *nBlock;

        m_blocks.push_back(std::move(nBlock));
        std::sort(m_blocks.begin(), m_blocks.end(), BlockedArray2<T, IDX>::block_cmp);

        return nBlockRef;
    }

    void merge(BlockedArray2 &&other) {
        for (uint i = 0; i < other.m_blocks.size(); ++i)
            m_blocks.push_back(std::move(other.m_blocks[i]));
        std::sort(m_blocks.begin(), m_blocks.end(), BlockedArray2<T, IDX>::block_cmp);
    }

//    const T &operator[](const IDX idx) const {
//        return at(idx, m_hint);
//    }
//
//    T &operator[](const IDX idx) {
//        return at(idx, m_hint);
//    }

    const T &unsafe_at(const IDX idx) const {
        return at(idx, m_hint);
    }

    T &unsafe_at(const IDX idx) {
        return at(idx, m_hint);
    }

    const T &at(const IDX idx, uint &hint) const {
        if (__builtin_expect(!_testBlock(idx, hint), false)) {
            hint = _findBlock(idx);
        }

        return m_blocks[hint]->at(idx);
    }

    T &at(const IDX idx, uint &hint) {
        if (__builtin_expect(!_testBlock(idx, hint), false)) {
            hint = _findBlock(idx);
        }

        return m_blocks[hint]->at(idx);
    }

    bool unsafe_contains(const IDX idx) const {
        return contains(idx, m_hint);
    }

    bool contains(const IDX idx, uint &hint) const {
        if (__builtin_expect(m_blocks[hint]->contains(idx), true))
            return true;

        uint block = _findBlock(idx, false);
        if (block != m_blocks.size() && m_blocks[block]->contains(idx)) {
            hint = block;
            return true;
        }

        return false;
    }

    IDX lowerBound() const {
        return m_blocks.front()->min();
    }

    IDX upperBound() const {
        return m_blocks.back()->max();
    }

public:
    typedef _detail::block_iterator<BlockedArray2, T> iterator;
    friend iterator;

    typedef _detail::block_iterator<const BlockedArray2, T> const_iterator;
    friend const_iterator;

    typedef _detail::range_type<const BlockedArray2, const_iterator> const_range_type;

    const_iterator begin() const {
        return const_iterator(*this, 0, lowerBound());
    }

    const_iterator end() const {
        return const_iterator(*this, m_blocks.size() - 1, upperBound());
    }

    iterator begin() {
        return iterator(*this, 0, lowerBound());
    }

    iterator end() {
        return iterator(*this, m_blocks.size() - 1, upperBound());
    }

    const_range_type range() const {
        return const_range_type(*this);
    }

protected:

    bool _testBlock(const IDX idx, const uint block) const {
        return m_blocks[block]->min() <= idx && idx < m_blocks[block]->max();
    }

    uint _findBlock(const IDX idx, bool doThrow = true) const {
        uint start = 0;
        uint end = m_blocks.size();
        uint it, count, step;
        count = end - start;

        while (count > 0) {
            step = count / 2;
            it = start + step;
            if (_testBlock(idx, it))
                return it;
            else {
                if (m_blocks[it]->min() < idx) {
                    start = ++it;
                    count -= step + 1;
                } else count = step;
            }
        }

        if (start == end && doThrow) {
            throw std::out_of_range(
                    "index: " + std::to_string(idx) + " no block found");
        }

        return start;
    }

    static bool block_cmp(const block_ptr &a, const block_ptr &b) {
        return a->min() < b->min();
    }

private:
    friend class boost::serialization::access;

    template<class Archive>
    void save(Archive &ar, __attribute__((unused)) const unsigned int version) const {
        // invoke serialization of the base class
        std::size_t n = m_blocks.size();
        ar << n;

        for (const auto & b : m_blocks) {
            ar << *b;
        }
    }

    template<class Archive>
    void load(Archive &ar, __attribute__((unused)) const unsigned int version) {
        // invoke serialization of the base class
        std::size_t n;
        ar >> n;
        m_blocks.resize(n);

        for (std::size_t i = 0; i < n; ++i) {
            ar >> *m_blocks[i];
        }
    }

    BOOST_SERIALIZATION_SPLIT_MEMBER()

protected:
    std::vector<block_ptr> m_blocks;
    mutable uint m_hint = 0;

};

template<class BA, typename IDX = std::size_t>
class BlockedArray2Handle {

public:
    BlockedArray2Handle(BA &container) : m_container(container) { }

    auto &operator[](const IDX idx) const {
        return m_container.at(idx, m_hint);
    }

    auto &operator[](const IDX idx) {
        return m_container.at(idx, m_hint);
    }

    const auto &at(const IDX idx) const {
        return m_container.at(idx, m_hint);
    }

    auto &at(const IDX idx) {
        return m_container.at(idx, m_hint);
    }

    bool contains(const IDX idx) const {
        return m_container.contains(idx, m_hint);
    }

private:
    BA &m_container;
    mutable uint m_hint = 0;
};