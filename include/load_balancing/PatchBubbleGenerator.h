#pragma once

#include <cassert>
#include <kdtree++/kdtree.hpp>
#include "utils/Generator.h"
#include "VectorOperations.h"

namespace LoadBalancing
{
template<uint D, typename Precision>
class PatchBubbleGenerator : public PointGenerator<D, Precision> {
public:
    PatchBubbleGenerator(size_t numBubbles)
     : mNumBubbles(numBubbles)
    {
        assert(numBubbles > 1);
    }
    
private:
    const size_t mNumBubbles;

protected:

    void _gen(dPoints<D, Precision> &points,
              const tIdType n,
              const dBox<D, Precision> &bounds,
              tGenerator &gen) const {
                  
        auto bubbleCenters = generateBubbleCenters(bounds, gen);
        auto bubbleRadiuses = calculateBubbleRadiuses(bubbleCenters.begin(), bubbleCenters.end());
                  
        std::normal_distribution<Precision> pointDist(0.0, 0.5);
        
        for(size_t i = 0; i < n; ++i) {
            size_t centerIndex = i * ((float)mNumBubbles / n);
            const auto& center = bubbleCenters[centerIndex];
            const auto& radius = bubbleRadiuses[centerIndex];
                  
            dVector<D, Precision> coords;
            do {
                std::generate(coords.begin(), coords.end(), [&gen, &pointDist]() -> auto {
                    return pointDist(gen);
                });
                coords = coords * radius + center;
            } while(!bounds.contains(coords));            
            points[i] = coords;
        }
    }
    
    std::vector<dVector<D, Precision>> generateBubbleCenters(const dBox<D, Precision> &bounds, tGenerator &gen) const {
        std::vector<dVector<D, Precision>> bubbleCenters(mNumBubbles);
        std::uniform_real_distribution<Precision> bubbleDist(0, 1);
        std::generate(bubbleCenters.begin(), bubbleCenters.end(), [&bounds, &gen, &bubbleDist]() -> auto
        {
            dVector<D, Precision> v;
            std::transform(bounds.low.begin(), bounds.low.end(), bounds.high.begin(), v.begin(),
                           [&gen, &bubbleDist](auto low, auto high) -> auto {
                               return low + (high - low) * bubbleDist(gen);
                        });
            return v;
        });
        return bubbleCenters;
    }
    
    template <typename It>
    std::vector<Precision> calculateBubbleRadiuses(It first, It last) const {
        std::vector<Precision> result(std::distance(first, last));
        std::transform(first, last, result.begin(),
                [first, last](const auto& center) {
                    std::vector<Precision> squaredDistances(std::distance(first, last));
                    std::transform(first, last, squaredDistances.begin(), [&center](const auto& other) {
                        return lenSquared(other - center);
                        
                    });
                    std::nth_element(squaredDistances.begin(), squaredDistances.begin() + 1, squaredDistances.end());
                    return std::sqrt(squaredDistances[1]) / 2;
                });
        return result;
    }
};
}
