//
// Created by dfunke on 5/7/15.
//

#pragma once

#include <type_traits>

template <typename T, typename std::enable_if<std::is_unsigned<T>::value>::type* = nullptr>
T nextPow2(const T x){
    T v = x-1;

    for(uint s = 1; s < sizeof(T)*8; s <<= 1){
        v |= v >> s;
    }

    return ++v;
}