#pragma once

#include <vector>
#include <algorithm>

template <class T>
class IdVector : public std::vector<T> {

public:
	T & operator[](uint idx){
		if(idx < std::vector<T>::size())
			return std::vector<T>::operator [](idx);
		else
			return *std::find(std::vector<T>::rbegin(), std::vector<T>::rend(), idx);
	}

	const T & operator[](uint idx) const {
		if(idx < std::vector<T>::size())
			return std::vector<T>::operator [](idx);
		else
			return *std::find(std::vector<T>::rbegin(), std::vector<T>::rend(), idx);
	}

};
