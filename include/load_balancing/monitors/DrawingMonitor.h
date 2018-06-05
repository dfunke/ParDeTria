#pragma once
#include <cstdlib>
#include <iostream>
#include <vector>
#include <tuple>
#include <functional>
#include <mutex>
#include "Painter.h"
#include "Monitor.h"

std::mutex mutex;

namespace LoadBalancing
{ 
	template <typename Precision>
    struct DrawingMonitor : public Monitor
    {
	    using PainterCreator = std::function<Painter<2, Precision>&(const std::string&)>;
	    
	    template <typename ForwardIt>
        DrawingMonitor(const dPoints<2, Precision>& points, ForwardIt colorBegin, ForwardIt  colorEnd,
                       std::mutex& mutex, PainterCreator createPainter)
          : mPoints(&points), mColors(colorBegin, colorEnd), mCreatePainter(std::move(createPainter)),
            mMutex(&mutex) {
        }
		
        void registerPartialTriangulation
	        (const std::vector<dSimplices<2, Precision>>& partialTriangulations,
	         const dSimplices<2, Precision>& borderTriangulation,
	         const dPoints<2, Precision>& points, const std::string& provenance)
	        {
				std::lock_guard<decltype(mutex)> lock(mutex);
		        auto& painter = mCreatePainter(provenance);
		        size_t c = 0;
		        for(const auto& triang : partialTriangulations) {
			        painter.setColor(mColors[c++ % mColors.size()]);
			        painter.draw(triang, points);
				}
				painter.setColor(tRGB(0,0,0));
			    painter.draw(borderTriangulation, points);
	        };
        
    private:
        const dPoints<2, Precision>* mPoints;
		std::vector<tRGB> mColors;
        PainterCreator mCreatePainter;
		std::mutex* mMutex;
    };
}
