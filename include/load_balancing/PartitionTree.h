#pragma once
#include <vector>
#include <variant>
#include <memory>
#include <algorithm>
#include "Geometry.h"
#include "../Partitioner.h"

namespace LoadBalancing
{    
    template <uint D, typename Precision>
    struct PartitionTree
    {
        using ChildContainer = std::vector<PartitionTree>;
        
        dBox<D, Precision> bounds;
        std::variant<ChildContainer, Point_Ids> attachment;
        //size_t numberOfPoints;
        void collect()
        {
            if(auto children = std::get_if<ChildContainer>(&attachment)) {
                Point_Ids newAttachment;
                for(auto& child : *children) {
                    child.collect();
                    auto& ids = std::get<Point_Ids>(child.attachment);
                    for(auto id : ids)
                        newAttachment.insert(id);
                }
                attachment = std::move(newAttachment);
            }
        }
        
        bool isLeaf() const {
            return std::holds_alternative<Point_Ids>(attachment);
        }
    };
        
    template <uint D, typename Precision>
    Partitioning<D, Precision> toOldPartitioning(typename PartitionTree<D, Precision>::ChildContainer children) {
        Partitioning<D, Precision> result;
        size_t id = 0;
        for(auto& childTree : children) {
            childTree.collect();
            auto ids = std::move(std::get<Point_Ids>(childTree.attachment));
            Partition<D, Precision> partition(ids.size());
            partition.points = std::move(ids);
            partition.bounds = std::move(childTree.bounds);
            partition.id = id++;
            result.push_back(std::move(partition));
        }
        return result;
    }
}
