#pragma once

#include <cassert>
#include "utils/Generator.h"
#include "VectorOperations.h"
#include "BoxUtils.h"

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
        auto bubbleRadiuses = calculateBubbleRadiuses(bubbleCenters.begin(), bubbleCenters.end(), bounds);
                  
        std::normal_distribution<Precision> coordsDist(0.0, 0.5);
        std::normal_distribution<Precision> pointsDist(1.0, 0.1);

        auto remainingPoints = n;
        size_t idx = 1;

        for(size_t c = 0; c < mNumBubbles; ++c){

            size_t pp = 0;

            if(mNumBubbles - c > 1)
                pp = (remainingPoints / (mNumBubbles - c)) * pointsDist(gen);
            else
                pp = remainingPoints;

            remainingPoints -= pp;

            for(size_t i = 0; i < pp; ++i) {
                const auto& center = bubbleCenters[c];
                const auto& radius = bubbleRadiuses[c];

                dVector<D, Precision> coords;
                do {
                    std::generate(coords.begin(), coords.end(), [&gen, &coordsDist]() -> auto {
                        return coordsDist(gen);
                    });
                    coords = coords * radius + center;
                } while(!bounds.contains(coords));
                points[idx++] = coords;
            }
        }

        ASSERT(remainingPoints == 0);
        ASSERT(idx == n + 1);
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
    std::vector<Precision> calculateBubbleRadiuses(It first, It last, const dBox<D, Precision>& boundary) const {
        std::vector<Precision> nnSquaredDistances(std::distance(first, last));
        std::transform(first, last, nnSquaredDistances.begin(),
                [first, last](const auto& center) {
                    std::vector<Precision> squaredDistances(std::distance(first, last));
                    std::transform(first, last, squaredDistances.begin(), [&center](const auto& other) {
                        return lenSquared(other - center);
                        
                    });
                    std::nth_element(squaredDistances.begin(), squaredDistances.begin() + 1, squaredDistances.end());
                    return squaredDistances[1];
                });
        
        std::vector<Precision> boundarySquaredDistances(std::distance(first, last));
        std::transform(first, last, boundarySquaredDistances.begin(),
                [first, last, &boundary](const auto& center) {
                    return lenSquared(vecToBoxBoundary<D, Precision>(center, boundary));
                });
        
        std::vector<Precision> result(std::distance(first, last));
        std::transform(nnSquaredDistances.begin(), nnSquaredDistances.end(), boundarySquaredDistances.begin(), result.begin(),
                       [](Precision nnSqDist, Precision boundSqDist) {
                           return std::sqrt(std::min(nnSqDist, boundSqDist)) / 2;
                    });
        return result;
    }
};
}
