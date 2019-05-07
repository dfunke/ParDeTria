#pragma once

#include "load_balancing/GridIntersectionChecker.h"
#include "load_balancing/GridIntersectionChecker.h"

const std::vector<tRGB> colors = {
    {1.0, 0.0, 0.0},
    {0.0, 1.0, 0.0},
    {0.0, 0.0, 1.0},
    {1.0, 1.0, 0.0},
    {0.0, 1.0, 1.0},
    {1.0, 0.0, 1.0},
};

template <typename Precision>
struct PartitionTreePainter
{
    PartitionTreePainter(Painter<2, Precision>& painter, const dPoints<2, Precision>& points)
     : painter(&painter), points(&points)
    {}
    
    void operator()(const lb::PartitionTree<2, Precision>& tree) {
        size_t partitions = 0;
        paintPartitionTree(tree, colors.cbegin(), partitions);
        
        /*painter->setColor(tRGB(0, 0, 0));
        for(const auto label : labels)
            painter->drawText(std::get<0>(label), std::get<1>(label));
        labels.clear();*/
    }
    

    std::vector<tRGB>::const_iterator paintPartitionTree(const lb::PartitionTree<2, Precision>& tree,
                                                   std::vector<tRGB>::const_iterator current,
                                                   size_t& currentNumPartitions) {
        if(current == colors.cend())
            current = colors.cbegin();
        
        if(tree.isLeaf()) {
            auto color = *current++;
            painter->setColor(color);
            dVector<2, Precision> center;
            size_t currentNumPoints = 0;
            for(auto id : std::get<Point_Ids>(tree.attachment)) {
                painter->draw((*points)[id]);
                center = center * ((Precision)currentNumPoints/(1 + currentNumPoints))
                    + (*points)[id].coords * ((Precision)1/(1 + currentNumPoints));
            }
            

			using GIC = lb::GridIntersectionChecker<2, Precision>;
			auto gridChecker = dynamic_cast<GIC*>(tree.intersectionChecker.get());
			if(!gridChecker) {
	            painter->drawBox(tree.intersectionChecker->bounds());
			} else {
				for(auto it = gridChecker->cellsBegin(); it != gridChecker->cellsEnd(); ++it) {
					auto hasNeighbour =
						[it, gridChecker] (const dIndex<2, Precision>& direction) -> bool {
							auto cell = *it;
							std::transform(cell.begin(), cell.end(), direction.begin(),
							               cell.begin(), [] (auto l, auto r) {
							               return l + r;
							               });
							return std::find(gridChecker->cellsBegin(), gridChecker->cellsEnd(),
							                 cell) != gridChecker->cellsEnd();
					};
					auto bounds = gridChecker->grid().boundsOf(*it);
					auto bottomLeft = bounds.low;
					auto topRight = bounds.high;
					auto bottomRight = bottomLeft;
					bottomRight[0] = topRight[0];
					auto topLeft = bottomLeft;
					topLeft[1] = topRight[1];
					if(hasNeighbour({0, 1})) {
						painter->drawLine(topLeft, topRight);
					}
					if(hasNeighbour({1, 0})) {
						painter->drawLine(topRight, bottomRight);
					}
					if(hasNeighbour({0, -1})) {
						painter->drawLine(bottomLeft, bottomRight);
					}
					if(hasNeighbour({-1, 0})) {
						painter->drawLine(topLeft, bottomLeft);
					}
				}
			}
            
            labels.push_back({std::to_string(currentNumPartitions++), center});
        } else {
            const auto& children = std::get<typename lb::PartitionTree<2, Precision>::ChildContainer>(tree.attachment);
            for(const auto& child : children) {
                current = paintPartitionTree(child, current, currentNumPartitions);
            }
        }
        return current;
    }
    
private:
    Painter<2, Precision>* painter;
    const dPoints<2, Precision>* points;
    std::vector<std::tuple<std::string, dVector<2, Precision>>> labels;
};
