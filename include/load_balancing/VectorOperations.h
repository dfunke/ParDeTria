#pragma once
#include <cmath>
#include <algorithm>
#include "Geometry.h"

template <size_t D, typename Precision, typename T>
dVector<D, Precision> operator*(const dVector<D, Precision>& left, T factor) {
    dVector<D, Precision> result;
    std::transform (left.cbegin(), left.cend(), result.begin(), [factor](Precision element) -> Precision {
        return element * factor;
    });
    return result;
}

template <size_t D, typename Precision, typename T>
dVector<D, Precision> operator*(T factor, const dVector<D, Precision>& left) {
    return left * factor;
}

template <size_t D, typename Precision>
dVector<D, Precision> operator+(const dVector<D, Precision>& left, const dVector<D, Precision>& right) {
    dVector<D, Precision> result;
    std::transform (left.cbegin(), left.cend(), right.cbegin(), result.begin(), std::plus<Precision>());
    return result;
}

template <size_t D, typename Precision>
dVector<D, Precision> operator-(const dVector<D, Precision>& vec) {
    return (-1) * vec;
}

template <size_t D, typename Precision>
dVector<D, Precision> operator-(const dVector<D, Precision>& left, const dVector<D, Precision>& right) {
    return left + (-right);
}

template <size_t D, typename Precision>
Precision scalarProduct(const dVector<D, Precision>& left, const dVector<D, Precision>& right) {
    return std::inner_product(left.cbegin(), left.cend(), right.cbegin(), 0.0);
}

template <size_t D, typename Precision>
Precision lenSquared(const dVector<D, Precision>& vec) {
    return scalarProduct(vec, vec);
}

template <size_t D, typename Precision>
Precision len(const dVector<D, Precision>& vec) {
    return std::sqrt(lenSquared(vec));
}

template <size_t D, typename Precision>
Precision distanceToPlane(const dVector<D, Precision>& point,
                          const dVector<D, Precision>& normal, Precision offset) {
			return scalarProduct(normal, point) - offset;
}
