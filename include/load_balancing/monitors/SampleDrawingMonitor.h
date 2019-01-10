#pragma once
#include <cstdlib>
#include <iostream>
#include <vector>
#include <tuple>
#include <functional>
#include <mutex>
#include "Painter.h"
#include "Monitor.h"


namespace LoadBalancing
{ 
	template <typename Precision>
    struct SampleDrawingMonitor : public Monitor
    {
	    using PainterCreator = std::function<Painter<2, Precision>&(const std::string&)>;
	    
	    template <typename ForwardIt>
        SampleDrawingMonitor(const dPoints<2, Precision>& points, ForwardIt colorBegin,
                             ForwardIt  colorEnd, std::mutex& mutex, PainterCreator createPainter)
          : mPoints(&points), mColors(colorBegin, colorEnd), mCreatePainter(std::move(createPainter)),
            mMutex(&mutex) {
        }
		
        void registerSampling(const Sampling<2, Precision>& sampling,
                              const std::string& provenance) {
	        std::lock_guard<decltype(*mMutex)> lock(*mMutex);
	        paintSampling(sampling, mCreatePainter(provenance));
		}
        
    private:
        const dPoints<2, Precision>* mPoints;
		std::vector<tRGB> mColors;
        PainterCreator mCreatePainter;
		std::mutex* mMutex;

		void paintSampling(const Sampling<2, Precision>& sampling,
						   Painter<2, Precision>& samplePainter)
		{
			const auto& graph = sampling.graph;
			for(size_t i = 0; i < graph.nodeRecords.size() - 1; ++i) {
				for(int k = graph.nodeRecords[i]; k < graph.nodeRecords[i + 1]; ++k) {
					size_t j = graph.adjacency[k];
					samplePainter.setColor(tRGB(0, 0, 0));
					samplePainter.drawLine(sampling.points[i], sampling.points[j]);
					
					auto center = 0.5 * (sampling.points[i] + sampling.points[j]);
					std::stringstream ss;
					ss << graph.edgeWeights[k];
					samplePainter.drawText(ss.str(), center, 10);
				}
				samplePainter.setColor(mColors[sampling.partition[i] % colors.size()]);
				samplePainter.draw(sampling.points[i]);
			}
		}
	};
}
