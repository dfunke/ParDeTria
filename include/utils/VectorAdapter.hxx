#pragma once

#include <vector>
#include <stdexcept>

template<typename T>
class VectorAdapter : public std::vector<T> {

public:
    typedef typename std::vector<T> vector;

public:

    VectorAdapter(const typename vector::size_type firstIdx = 0)
            : vector(), m_firstIdx(firstIdx) { }

public:
    bool contains(const typename vector::size_type &i) const {

        return (T::isFinite(i) ? _finIdx(i) < finite_size() : _infIdx(i)) < vector::size();
    }

    T &operator[](const typename vector::size_type i) {
        if (__builtin_expect(T::isFinite(i), true))
            return vector::operator[](_finIdx(i));
        else
#ifndef NDEBUG
            return vector::at(_infIdx(i));
#else
      return vector::operator[](_infIdx(i));
#endif
    }

    const T &operator[](const typename vector::size_type i) const {
        if (__builtin_expect(T::isFinite(i), true))
            return vector::operator[](_finIdx(i));
        else
#ifndef NDEBUG
            return vector::at(_infIdx(i));
#else
      return vector::operator[](_infIdx(i));
#endif
    }

    T &at(const typename vector::size_type i) {
        if (__builtin_expect(T::isFinite(i), true))
            return vector::at(_finIdx(i));
        else
            return vector::at(_infIdx(i));
    }

    const T &at(const typename vector::size_type i) const {
        if (__builtin_expect(T::isFinite(i), true))
            return vector::at(_finIdx(i));
        else
            return vector::at(_infIdx(i));
    }

    template<class Return, class Container>
    const Return filter(const Container &ids) const {
        return _filter<Return>(ids.begin(), ids.end(), ids.size());
    }

    template<class Return, class InputIt>
    const Return filter(const InputIt &first, const InputIt &last) const {
        return _filter<Return>(first, last, std::distance(first, last));
    }

    const typename vector::size_type finite_size() const {
        return vector::size() - T::nINF;
    }

private:
    inline typename vector::size_type
    _infIdx(const typename vector::size_type i) const {
        return vector::size() - T::nINF + T::infIndex(i);
    }

    inline typename vector::size_type
    _finIdx(const typename vector::size_type i) const {
        return i - m_firstIdx;
    }

    template<class Return, class InputIt>
    const Return _filter(const InputIt &first, const InputIt &last,
                         const uint n) const {
        Return res;
        res.reserve(n);

        for (auto it = first; it != last; ++it)
            res.emplace_back(at(*it));

        return res;
    }

private:
    typename vector::size_type m_firstIdx;
};
