#pragma once

#include <cassert>
#include "utils/Generator.h"
#include "VectorOperations.h"
#include "BoxUtils.h"

namespace LoadBalancing
{
template<uint D, typename Precision>
class MaliciousBubbleGenerator : public PointGenerator<D, Precision> {
public:
    MaliciousBubbleGenerator(size_t numBubbles)
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
                  
        auto bubbleCenters = generateBubbleCenters(bounds);
        auto bubbleRadiuses = calculateBubbleRadiuses(bubbleCenters.begin(), bubbleCenters.end(), bounds);
                  
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
    
    std::vector<dVector<D, Precision>> generateBubbleCenters(const dBox<D, Precision> &bounds) const {
        std::vector<dVector<D, Precision>> bubbleCenters;
        bubbleCenters.reserve(mNumBubbles);

        std::uniform_real_distribution<Precision> bubbleDist(0, 1);

        uint bubblesPerDim = static_cast<uint>(std::ceil(std::pow(mNumBubbles, 1.0 / D)));

        auto inc = bounds.high;
        for(uint d = 0; d < D; ++d) {
            inc[d] -= bounds.low[d];
            inc[d] /= bubblesPerDim;
        }

        //x = id % chunks_per_dim_;
        //y = (id / chunks_per_dim_) % chunks_per_dim_;
        //z = id / (chunks_per_dim_ * chunks_per_dim_);


        for(uint i = 0; i < mNumBubbles; ++i){
            auto center = bounds.low;

            for(uint d = 0; d < D; ++d){
                center[d] += inc[d] * ((i / static_cast<uint>(std::pow(bubblesPerDim, (d)))) % bubblesPerDim);
            }
            bubbleCenters.push_back(center);
        }

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
                           return std::sqrt(std::min(nnSqDist, boundSqDist)) / 4;
                    });
        return result;
    }
};
}
