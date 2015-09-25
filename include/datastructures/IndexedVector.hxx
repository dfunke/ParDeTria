#pragma once

#include <unordered_map>
#include <vector>
#include <iterator>

#include "utils/TBB_Containers.h"
#include "utils/Misc.h"
#include "utils/Timings.h"

#include <csignal>

template<typename V>
class IndexedVector : private tbb::concurrent_vector<V> {

public:
    typedef typename tbb::concurrent_vector<V> base;

public:
    V &operator[](typename base::size_type idx) {
        PROFILER_INC("IndexedVector_access");

        try {
            return base::at(idx);
        } catch (std::exception &e) {
            std::cerr << e.what() << std::endl;
            std::cerr << "size: " << base::size() << " index: " << idx << std::endl;
            raise(SIGINT);
            throw e;
        }
    }

    const V &operator[](typename base::size_type idx) const {
        PROFILER_INC("IndexedVector_access");

        try {
            return base::at(idx);
        } catch (std::exception &e) {
            std::cerr << e.what() << std::endl;
            std::cerr << "size: " << base::size() << " index: " << idx << std::endl;
            raise(SIGINT);
            throw e;
        }
    }

    V &at(typename base::size_type idx) {
        PROFILER_INC("IndexedVector_access");

        try {
            return base::at(idx);
        } catch (std::exception &e) {
            std::cerr << e.what() << std::endl;
            std::cerr << "size: " << base::size() << " index: " << idx << std::endl;
            raise(SIGINT);
            throw e;
        }
    }

    const V &at(typename base::size_type idx) const {
        PROFILER_INC("IndexedVector_access");

        try {
            return base::at(idx);
        } catch (std::exception &e) {
            std::cerr << e.what() << std::endl;
            std::cerr << "size: " << base::size() << " index: " << idx << std::endl;
            raise(SIGINT);
            throw e;
        }
    }

    void insert(const V &value) {
        PROFILER_INC("IndexedVector_insert");

        base::push_back(value);
    }

    template<class InputIt>
    void insert(const InputIt &first, const InputIt &last) {
        for (auto it = first; it != last; ++it) {
            insert(*it);
        }
    }

    void reserve(std::size_t s) {
        PROFILER_INC("IndexedVector_reserve");

        base::reserve(nextPow2(s));
    }

    void grow(std::size_t s) {
        PROFILER_INC("IndexedVector_reserve");

        base::grow_to_at_least(s);
    }

    template<class Container>
    const IndexedVector project(const Container &ids) const {
        return project(ids.begin(), ids.end());
    }

    template<class InputIt>
    const IndexedVector project(const InputIt &first, const InputIt &last) const {
        IndexedVector res;

        for (auto it = first; it != last; ++it)
            res.insert(operator[](*it));

        return res;
    }

    template<class Return, class Container>
    const Return project(const Container &ids) const {
        return project<Return>(ids.begin(), ids.end());
    }

    template<class Return, class InputIt>
    const Return project(const InputIt &first, const InputIt &last) const {
        Return res;

        for (auto it = first; it != last; ++it)
            res.insert(operator[](*it));

        return res;
    }

    auto size() const {
        return base::size();
    }

    auto begin() {
        PROFILER_INC("IndexedVector_begin");

        return base::begin();
    }

    auto end() {
        return base::end();
    }

    auto begin() const {
        PROFILER_INC("IndexedVector_begin");

        return base::begin();
    }

    auto end() const {
        return base::end();
    }
};
