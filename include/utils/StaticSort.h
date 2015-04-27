#pragma once

#include <array>

template<typename T, std::size_t N>
void static_insertion_sort(std::array<T, N> &arr) {
    for (unsigned short i = 1; i < N; ++i) {
        double value = arr[i];
        short hole = i;

        for (; hole > 0 && value < arr[hole - 1]; --hole)
            arr[hole] = arr[hole - 1];

        arr[hole] = value;
    }
}
