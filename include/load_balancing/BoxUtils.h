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
    
    template <uint D, typename Precision>
    void enlargeBoxAroundVector(dBox<D, Precision>& box, const dVector<D, Precision>& v) {
        auto overlap = [](const dVector<D, Precision>& left, const dVector<D, Precision>& right) -> dVector<D, Precision> {
            dVector<D, Precision> result;
            for(size_t d = 0; d < D; ++d) {
                result[d] = std::max(left[d], right[d]);
            }
            return result;
        };
        
        auto underlap = [](const dVector<D, Precision>& left, const dVector<D, Precision>& right) -> dVector<D, Precision> {
            dVector<D, Precision> result;
            for(size_t d = 0; d < D; ++d) {
                result[d] = std::min(left[d], right[d]);
            }
            return result;
        };
        
        for(size_t d = 0; d < D; ++d) {
            box.low = underlap(box.low, v);
            box.high = overlap(box.high, v);
        }
    }
}
