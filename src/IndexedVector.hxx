#pragma once

#include <map>
#include <vector>
#include <iterator>

template <typename V, typename K = uint>
class IndexedVector: public std::map<K, V> {

public:

	typedef typename std::map<K, V> map;

public:

	struct const_iterator; //forward declare for friend declaration

	struct iterator: public std::iterator<std::bidirectional_iterator_tag, V> {

		friend struct const_iterator;

	public:

		iterator() {
		}
		iterator(typename map::iterator j) :
				i(j) {
		}

		iterator& operator++() {
			++i;
			return *this;
		}

		iterator operator++(int) {
			auto tmp = *this;
			++(*this);
			return tmp;
		}

		iterator& operator--() {
			--i;
			return *this;
		}

		iterator operator--(int) {
			auto tmp = *this;
			--(*this);
			return tmp;
		}

		iterator & operator=(const iterator & other) {
			i = other.i;

			return *this;
		}

		bool operator==(iterator j) const {
			return i == j.i;
		}
		bool operator!=(iterator j) const {
			return !(*this == j);
		}

		V & operator*() {
			return i->second;
		}

		V * operator->() {
			return &i->second;
		}

	protected:
		typename map::iterator i;
	};

	struct const_iterator: public std::iterator<std::bidirectional_iterator_tag, V> {

	public:

		const_iterator() {
		}

		const_iterator(typename map::const_iterator j) :
				i(j) {
		}

		const_iterator(iterator & o) :
				i(o.i) {
		}

		const_iterator& operator++() {
			++i;
			return *this;
		}

		const_iterator operator++(int) {
			auto tmp = *this;
			++(*this);
			return tmp;
		}

		const_iterator& operator--() {
			--i;
			return *this;
		}

		const_iterator operator--(int) {
			auto tmp = *this;
			--(*this);
			return tmp;
		}

		const_iterator & operator=(const const_iterator & other) {
			i = other.i;

			return *this;
		}

		const_iterator & operator=(const iterator & other) {
			i = other.i;

			return *this;
		}

		bool operator==(const_iterator j) const {
			return i == j.i;
		}

		bool operator!=(const_iterator j) const {
			return !(*this == j);
		}

		const V & operator*() const {
			return i->second;
		}

		const V * operator->() const {
			return &i->second;
		}

	protected:
		typename map::const_iterator i;
	};

public:
	V & operator[](uint idx) {
		return map::operator [](idx);
	}

	const V & operator[](uint idx) const {
		return map::at(idx);
	}

	void push_back(const V& value) {
		map::insert(std::make_pair(value.id, value));
	}

	template<class InputIt>
	void insert(InputIt first, InputIt last) {
		for(; first != last; ++first){
			push_back(*first);
		}
	}

	bool contains(const V & value) const {
		return map::count(value.id) == 1;
	}

	bool contains(const K & key) const {
		return map::count(key) == 1;
	}

	IndexedVector project(const std::vector<K> & ids) const {
		IndexedVector res;

		for(const auto & id : ids)
			res.push_back(operator[](id));

		return res;
	}

	iterator begin() { return iterator(map::begin()); }

	iterator end() { return iterator(map::end()); }

	const_iterator begin() const { return const_iterator(map::begin());	}

	const_iterator end() const { return const_iterator(map::end());	}
};