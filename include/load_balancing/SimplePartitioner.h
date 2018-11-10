#pragma once

#include <cassert>
#include <memory>
#include "Partitioner.h"
#include "BoxUtils.h"
#include "load_balancing/BoundsIntersectionChecker.h"

namespace LoadBalancing
{
    template <uint D, typename Precision, typename MonitorT>
    struct SimplePartitioner : public Partitioner<D, Precision, MonitorT>
    {        
        PartitionTree<D, Precision> partition(const dBox<D, Precision>& bounds,
                                        dPoints<D, Precision>& points,
                                        const Point_Ids& pointIds) override
        {
            auto x = (bounds.low[0] + bounds.high[0])/2;            
            auto boxPair = splitBox(bounds, 0, x);
            dBox<D, Precision> leftBounds = std::move(std::get<0>(boxPair));
            dBox<D, Precision> rightBounds = std::move(std::get<1>(boxPair));
            
            Point_Ids leftIds, rightIds;
            for(auto id : pointIds) {
                const auto& point = points[id];
                
                if(dPoint<D, Precision>::isFinite(id)) {
                    if(leftBounds.contains(point.coords))
                        leftIds.insert(id);
                    else
                        rightIds.insert(id);
                }
            }
        
            for (tIdType k = dPoint<D, Precision>::cINF; k != 0; ++k) {
                leftIds.insert(k);
                rightIds.insert(k);
            }
            
            PartitionTree<D, Precision> leftSubtree;
            leftSubtree.intersectionChecker = std::make_unique<BoundsIntersectionChecker<D, Precision>>(leftBounds);
            leftSubtree.attachment = std::move(leftIds);

            PartitionTree<D, Precision> rightSubtree;
            rightSubtree.intersectionChecker = std::make_unique<BoundsIntersectionChecker<D, Precision>>(rightBounds);
            rightSubtree.attachment = std::move(rightIds);
            
            PartitionTree<D, Precision> tree;
            tree.intersectionChecker = std::make_unique<BoundsIntersectionChecker<D, Precision>>(bounds);
            typename PartitionTree<D, Precision>::ChildContainer children;
            children.push_back(std::move(leftSubtree));
            children.push_back(std::move(rightSubtree));
            tree.attachment = std::move(children);
            
            
            return tree;
        }
        
        std::string info() const override
        {
            return "simple partitioner";
        }
    };
}
