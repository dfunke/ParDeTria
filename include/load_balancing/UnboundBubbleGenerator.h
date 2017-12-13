#pragma once

#include "utils/Generator.h"

namespace LoadBalancing
{
template<uint D, typename Precision>
class UnboundBubbleGenerator : public PointGenerator<D, Precision> {
public:
    UnboundBubbleGenerator(size_t numBubbles, Precision bubbleRadius)
     : mNumBubbles(numBubbles), mBubbleRadius(bubbleRadius)
    {}
    
private:
    size_t mNumBubbles;
    Precision mBubbleRadius;

protected:

    void _gen(dPoints<D, Precision> &points,
              const tIdType n,
              const dBox<D, Precision> &bounds,
              tGenerator &gen) const {


        std::uniform_real_distribution<Precision> bubbleDist(0, 1);
        std::uniform_real_distribution<Precision> pointDist(-mBubbleRadius, mBubbleRadius);
        
        std::vector<dVector<D, Precision>> bubbleCenters(mNumBubbles);
        std::generate(bubbleCenters.begin(), bubbleCenters.end(), [&bounds, &bubbleDist, &gen]() -> auto
        {
            dVector<D, Precision> v;
            std::transform(bounds.low.begin(), bounds.low.end(), bounds.high.begin(), v.begin(), [&bubbleDist, &gen](auto low, auto high) -> auto {
                return low + (high - low) * bubbleDist(gen);
            });
            return v;
        });
        
        for(size_t i = 0; i < n; ++i) {
            size_t centerIndex = i * ((float)mNumBubbles / n);
            
            dVector<D, Precision> coords;
            std::generate(coords.begin(), coords.end(), [&pointDist, &gen]() -> auto {
                return pointDist(gen);
            });
            
            points[i] = bubbleCenters[centerIndex] + coords;
        }
    }
};
}
