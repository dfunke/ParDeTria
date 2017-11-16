#pragma once

#include <cassert>
#include "Geometry.h"
#include "load_balancing/VectorOperations.h"

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
    
    template <typename It>
    dBox<std::tuple_size<typename It::value_type>::value, typename It::value_type::value_type> makeBoundingBox(It vbegin, It vend) {
        constexpr uint D = std::tuple_size<typename It::value_type>::value;
        using Precision = typename It::value_type::value_type;
        
        dBox<D, Precision> result;
        if(vbegin != vend){
            result.low = *vbegin++;
            result.high = result.low;
            while(vbegin != vend){
                enlargeBoxAroundVector<D, Precision>(result, *vbegin++);
            }
        }
        return result;
    }
    
    template <uint D, typename Precision>
    dVector<D, Precision> vecToLowerQuadrant(dVector<D, Precision> v) {        
        std::transform(begin(v), end(v), begin(v), [](Precision x) {
            return std::max(x, 0.0);
        });
        return v;
    }
    
    template <uint D, typename Precision>
    dVector<D, Precision> vecToCenteredBox(dVector<D, Precision> v, const dVector<D, Precision>& highCorner) {
        std::transform(begin(v), end(v), begin(v), [](Precision p) {
            return std::abs(p);
        });
        return vecToLowerQuadrant<D, Precision>(v - highCorner);
    }
    
    template <uint D, typename Precision>
    dVector<D, Precision> vecToBox(const dVector<D, Precision>& v, const dBox<D, Precision>& box) {
        auto center = (box.high + box.low) * Precision(0.5);
        return vecToCenteredBox<D, Precision>(v - center, box.high - center);
    }
}
