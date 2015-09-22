//
// Created by dfunke on 5/7/15.
//

#pragma once

#include <type_traits>
#include <limits.h>

template<typename T, typename std::enable_if<std::is_unsigned<T>::value>::type * = nullptr>

T nextPow2(const T x) {
    T v = x - 1;

    for (uint s = 1; s < sizeof(T) * CHAR_BIT; s <<= 1) {
        v |= v >> s;
    }

    return ++v;
}

template<typename T, typename std::enable_if<std::is_unsigned<T>::value>::type * = nullptr>

T log2(const T _x) {
    T x = _x;
    T y = 0;
    while (x >>= 1) ++y;
    //asm ( "\tbsr %1, %0\n" : "=r"(y) : "r" (x));
    return y;
}

template<typename T, typename std::enable_if<std::is_unsigned<T>::value>::type * = nullptr>

T popcount(const T _x) {
    T v = _x;
    v = v - ((v >> 1) & (T)
    ~(T) 0 / 3);                           // temp
    v = (v & (T)
    ~(T) 0 / 15 * 3) +((v >> 2) & (T)
    ~(T) 0 / 15 * 3);      // temp
    v = (v + (v >> 4)) & (T)
    ~(T) 0 / 255 * 15;                      // temp
    return (T)(v * ((T)
    ~(T) 0 / 255)) >> (sizeof(T) - 1) * CHAR_BIT; // count
}