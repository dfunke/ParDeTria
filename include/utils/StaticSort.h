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

template<typename T, typename D, std::size_t N>
void static_insertion_sort_with_data(std::array<T, N> &arr, std::array<D, N> &dat) {
    static_assert(N < std::numeric_limits<unsigned short>::max(), "StaticSort: array to large");

    for (unsigned short i = 1; i < N; ++i) {
        T value = arr[i];
        D data = dat[i];
        short hole = i;

        for (; hole > 0 && value < arr[hole - 1]; --hole) {
            arr[hole] = arr[hole - 1];
            dat[hole] = dat[hole - 1];
        }

        arr[hole] = value;
        dat[hole] = data;
    }
}