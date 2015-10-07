#pragma once

#include <vector>
#include <stdexcept>

#include <tbb/blocked_range.h>

#include "utils/Timings.h"

template<typename T>
class VectorAdapter2 : private std::vector<T> {

public:
    typedef typename std::vector<T> vector;

public:

    VectorAdapter2(const typename vector::size_type offset = 0)
            : vector(), m_offset(offset) { }

public:
    bool contains(const typename vector::size_type &i) const {
        PROFILER_INC("VectorAdapter2_contains");

        return (__builtin_expect(T::isFinite(i), true))
               ? (i >= m_offset && _finIdx(i) < finite_size())
               : (_infIdx(i) < T::nINF);
    }

    T &operator[](const typename vector::size_type i) {
        PROFILER_INC("VectorAdapter2_access");

#ifndef NDEBUG
        return at(i);
#else
        if (__builtin_expect(T::isFinite(i), true))
            return vector::operator[](_finIdx(i));
        else
      return m_infinites[_infIdx(i)];
#endif
    }

    const T &operator[](const typename vector::size_type i) const {
        PROFILER_INC("VectorAdapter2_access");

#ifndef NDEBUG
        return at(i);
#else
        if (__builtin_expect(T::isFinite(i), true))
            return vector::operator[](_finIdx(i));
        else
      return m_infinites[_infIdx(i)];
#endif
    }

    T &at(const typename vector::size_type i) {
        PROFILER_INC("VectorAdapter2_access");

        if (__builtin_expect(T::isFinite(i), true)) {
            ASSERT(m_offset <= i && _finIdx(i) < finite_size());
            return vector::at(_finIdx(i));
        } else {
            ASSERT(_infIdx(i) < T::nINF);
            return m_infinites.at(_infIdx(i));
    }
    }

    const T &at(const typename vector::size_type i) const {
        PROFILER_INC("VectorAdapter2_access");

        if (__builtin_expect(T::isFinite(i), true)) {
            ASSERT(m_offset <= i && _finIdx(i) < finite_size());
            return vector::at(_finIdx(i));
        } else {
            ASSERT(_infIdx(i) < T::nINF);
            return m_infinites.at(_infIdx(i));
    }
    }

    template<class Return, class Container>
    const Return filter(const Container &ids) const {
        return _filter<Return>(ids.begin(), ids.end(), ids.size());
    }

    template<class Return, class InputIt>
    const Return filter(const InputIt &first, const InputIt &last) const {
        return _filter<Return>(first, last, std::distance(first, last));
    }

    const typename vector::size_type size() const {
        return vector::size() + T::nINF;
    }

    const typename vector::size_type finite_size() const {
        return vector::size();
    }

    typename vector::iterator begin() {
        PROFILER_INC("VectorAdapter2_begin");

        return vector::begin();
    }

    typename vector::iterator end() {
        return vector::end();
    }

    typename vector::const_iterator begin() const {
        PROFILER_INC("VectorAdapter2_begin");

        return vector::begin();
    }

    typename vector::const_iterator end() const {
        return vector::end();
    }

//    template<typename... _Args>
//    void emplace_back(_Args &&... args) {
//        vector::emplace_back(args...);
//    }

    bool reserveUpToIdx(const typename vector::size_type &size) {
        if(vector::size() == 0 || size - m_offset > vector::size()) {
            vector::resize(size - m_offset);
            return true;
        }

        return false;
    }

    void reserveEntries(const typename vector::size_type &size) {
        if (size > vector::size())
            vector::resize(size);
    }

    void offset(const typename vector::size_type &offset) {
        m_offset = offset;
    }

    std::size_t offset() const {
        return m_offset;
    }

    template<class InputIt, class PosIt>
    void insert(PosIt pos, InputIt begin, InputIt end) {
        vector::insert(pos, begin, end);
    }

    auto range() const {
        return tbb::blocked_range<typename vector::size_type>(m_offset, m_offset + finite_size());
    }

private:
    inline typename vector::size_type
    _infIdx(const typename vector::size_type i) const {
        return T::infIndex(i);
    }

    inline typename vector::size_type
    _finIdx(const typename vector::size_type i) const {
        return i - m_offset;
    }

    template<class Return, class InputIt>
    const Return _filter(const InputIt &first, const InputIt &last,
                         const typename vector::size_type n) const {
        Return res;
        res.reserve(n);

        for (auto it = first; it != last; ++it)
            res.emplace_back(at(*it));

        return res;
    }

private:
    friend class boost::serialization::access;

    template<class Archive>
    void serialize(Archive &ar, __attribute((unused)) const unsigned int version) {
        ar &boost::serialization::base_object<vector>(*this);
        ar &m_offset;
        ar &m_infinites;
    }

private:
    typename vector::size_type m_offset;
    std::array<T, T::nINF> m_infinites;
};
