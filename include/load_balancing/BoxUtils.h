#pragma once

#include <cassert>
#include "Geometry.h"

namespace LoadBalancing
{
    template <uint D, typename Precision>
    std::tuple<dBox<D, Precision>, dBox<D, Precision>> splitBox(const dBox<D, Precision>& bounds, uint dimension, Precision splitPosition) {
        assert(dimension < D);
        
        auto leftHigh = bounds.high;
        leftHigh[dimension] = splitPosition;
        auto rightLow = bounds.low;
        rightLow[dimension] = splitPosition;
        
        dBox<D, Precision> leftBounds(bounds.low, leftHigh);
        dBox<D, Precision> rightBounds(rightLow, bounds.high);
        
        return std::make_tuple(std::move(leftBounds), std::move(rightBounds));
    }
}
