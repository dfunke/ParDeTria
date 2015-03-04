#pragma once

#include <unordered_map>
#include <vector>
#include <iterator>

template <typename V, typename K = uint>
class IndexedVector : public std::unordered_map<K, V> {

public:
  typedef typename std::unordered_map<K, V> map;

public:
  template <class IT>
  struct const_iterator; // forward declare for friend declaration

  template <class IT>
  struct iterator : public std::iterator<std::bidirectional_iterator_tag, V> {

    friend struct const_iterator<IT>;

  public:
    iterator() {}
    iterator(IT j) : i(j) {}

    iterator &operator++() {
      ++i;
      return *this;
    }

    iterator operator++(int) {
      auto tmp = *this;
      ++(*this);
      return tmp;
    }

    iterator &operator--() {
      --i;
      return *this;
    }

    iterator operator--(int) {
      auto tmp = *this;
      --(*this);
      return tmp;
    }

    iterator &operator=(const iterator &other) {
      i = other.i;

      return *this;
    }

    bool operator==(iterator j) const { return i == j.i; }
    bool operator!=(iterator j) const { return !(*this == j); }

    V &operator*() { return i->second; }

    V *operator->() { return &i->second; }

  protected:
    IT i;
  };

  template <class IT>
  struct const_iterator
      : public std::iterator<std::bidirectional_iterator_tag, V> {

  public:
    const_iterator() {}

    const_iterator(IT j) : i(j) {}

    const_iterator(iterator<IT> &o) : i(o.i) {}

    const_iterator &operator++() {
      ++i;
      return *this;
    }

    const_iterator operator++(int) {
      auto tmp = *this;
      ++(*this);
      return tmp;
    }

    const_iterator &operator--() {
      --i;
      return *this;
    }

    const_iterator operator--(int) {
      auto tmp = *this;
      --(*this);
      return tmp;
    }

    const_iterator &operator=(const const_iterator &other) {
      i = other.i;

      return *this;
    }

    const_iterator &operator=(const iterator<IT> &other) {
      i = other.i;

      return *this;
    }

    bool operator==(const_iterator j) const { return i == j.i; }

    bool operator!=(const_iterator j) const { return !(*this == j); }

    const V &operator*() const { return i->second; }

    const V *operator->() const { return &i->second; }

  protected:
    IT i;
  };

  template <class IT>
  struct const_key_iterator
      : public std::iterator<std::bidirectional_iterator_tag, K> {

  public:
    const_key_iterator() {}

    const_key_iterator(IT j) : i(j) {}

    const_key_iterator &operator++() {
      ++i;
      return *this;
    }

    const_key_iterator operator++(int) {
      auto tmp = *this;
      ++(*this);
      return tmp;
    }

    const_key_iterator &operator--() {
      --i;
      return *this;
    }

    const_key_iterator operator--(int) {
      auto tmp = *this;
      --(*this);
      return tmp;
    }

    bool operator==(const_key_iterator j) const { return i == j.i; }

    bool operator!=(const_key_iterator j) const { return !(*this == j); }

    const K &operator*() const { return i->first; }

    const K *operator->() const { return &i->first; }

  protected:
    IT i;
  };

public:
  V &operator[](uint idx) { return map::operator[](idx); }

  const V &operator[](uint idx) const { return map::at(idx); }

  void insert(const V &value) { map::insert(std::make_pair(value.id, value)); }

  template <class InputIt>
  void insert(const InputIt &first, const InputIt &last) {
    for (auto it = first; it != last; ++it) {
      insert(*it);
    }
  }

  bool contains(const V &value) const { return map::count(value.id) == 1; }

  bool contains(const K &key) const { return map::count(key) == 1; }

  template <class Container>
  const IndexedVector project(const Container &ids) const {
    return project(ids.begin(), ids.end());
  }

  template <class InputIt>
  const IndexedVector project(const InputIt &first, const InputIt &last) const {
    IndexedVector res;

    for (auto it = first; it != last; ++it)
      res.insert(operator[](*it));

    return res;
  }

  template <class Return, class Container>
  const Return project(const Container &ids) const {
    return project<Return>(ids.begin(), ids.end());
  }

  template <class Return, class InputIt>
  const Return project(const InputIt &first, const InputIt &last) const {
    Return res;

    for (auto it = first; it != last; ++it)
      res.insert(operator[](*it));

    return res;
  }

  iterator<typename map::iterator> begin() {
    return iterator<typename map::iterator>(map::begin());
  }
  iterator<typename map::iterator> end() {
    return iterator<typename map::iterator>(map::end());
  }

  const_iterator<typename map::const_iterator> begin() const {
    return const_iterator<typename map::const_iterator>(map::begin());
  }
  const_iterator<typename map::const_iterator> end() const {
    return const_iterator<typename map::const_iterator>(map::end());
  }

  const_key_iterator<typename map::const_iterator> begin_keys() const {
    return const_key_iterator<typename map::const_iterator>(map::begin());
  }
  const_key_iterator<typename map::const_iterator> end_keys() const {
    return const_key_iterator<typename map::const_iterator>(map::end());
  }

  iterator<typename map::local_iterator> begin(const uint i) {
    return iterator<typename map::local_iterator>(map::begin(i));
  }
  iterator<typename map::local_iterator> end(const uint i) {
    return iterator<typename map::local_iterator>(map::end(i));
  }
};
