#pragma once

#include <array>
#include <limits>

template<typename T, std::size_t N>
void static_insertion_sort(std::array<T, N> &arr) {
    static_assert(N < std::numeric_limits<unsigned short>::max(), "StaticSort: array to large");

    for (unsigned short i = 1; i < N; ++i) {
        T value = arr[i];
        short hole = i;

        for (; hole > 0 && value < arr[hole - 1]; --hole)
            arr[hole] = arr[hole - 1];

        arr[hole] = value;
    }
}
