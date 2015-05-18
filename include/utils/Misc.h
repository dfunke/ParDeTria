//
// Created by dfunke on 5/7/15.
//

#pragma once

#include <type_traits>
#include <limits.h>

template <typename T, typename std::enable_if<std::is_unsigned<T>::value>::type* = nullptr>
T nextPow2(const T x){
    T v = x-1;

    for(uint s = 1; s < sizeof(T)*CHAR_BIT; s <<= 1){
        v |= v >> s;
    }

    return ++v;
}

template<typename T, typename std::enable_if<std::is_unsigned<T>::value>::type* = nullptr>
T log2(const T _x) {
    T x = _x;
    T y = 0;
    while (x >>= 1) ++y;
    //asm ( "\tbsr %1, %0\n" : "=r"(y) : "r" (x));
    return y;
}