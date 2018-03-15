#pragma once

#include <memory>
#include "../Partitioner.h"
#include "load_balancing/Partitioner.h"
#include "load_balancing/BoundsIntersectionChecker.h"

namespace LoadBalancing
{
    template <uint D, typename Precision>
    struct OldPartitionerPartitioner : public Partitioner<D, Precision>
    {
        OldPartitionerPartitioner(std::unique_ptr<::Partitioner<D, Precision>> splitter, size_t maxRecursions, size_t baseCutoff)
            : splitter(std::move(splitter)), maxRecursions(maxRecursions), baseCutoff(baseCutoff) { }
        
        PartitionTree<D, Precision> partition(const dBox<D, Precision>& bounds,
                                        const dPoints<D, Precision>& points,
                                        const Point_Ids& pointIds) override
        {
            return buildTree(bounds, points, pointIds, maxRecursions);
        }
        
        std::string info() const override
        {
            return "old partitioner(" + splitter->to_string() + ")";
        }
        
    private:
        PartitionTree<D, Precision> buildTree(const dBox<D, Precision>& bounds,
                                        const dPoints<D, Precision>& points,
                                        const Point_Ids& pointIds,
                                        size_t remainingRecursions) const
        {
            Partitioning<D, Precision> partitioning = splitter->partition(pointIds, points, "");
            
            typename PartitionTree<D, Precision>::ChildContainer children;
            for(auto& partition : partitioning) {
                PartitionTree<D, Precision> subtree;
                if(remainingRecursions > 0 && pointIds.size() > baseCutoff) {
                    subtree = buildTree(partition.bounds, points,
                                        partition.points, remainingRecursions - 1);
                    
                } else {
                    subtree.attachment = std::move(partition.points);
            		subtree.intersectionChecker = std::make_unique<BoundsIntersectionChecker<D, Precision>>(partition.bounds);
                }
                children.push_back(std::move(subtree));
            }
            
            PartitionTree<D, Precision> result;
            result.attachment = std::move(children);
            result.intersectionChecker = std::make_unique<BoundsIntersectionChecker<D, Precision>>(bounds);
            return result;
        }
        
        std::unique_ptr<::Partitioner<D, Precision>> splitter;
        size_t maxRecursions;
        size_t baseCutoff;
    };
}
