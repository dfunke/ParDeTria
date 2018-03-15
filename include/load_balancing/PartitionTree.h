#pragma once
#include <vector>
#include <variant>
#include <memory>
#include <algorithm>
#include <tuple>
#include "Geometry.h"
#include "../Partitioner.h"
#include "load_balancing/IntersectionChecker.h"

namespace LoadBalancing
{
    template <uint D, typename Precision>
    struct PartitionTree
    {
        using ChildContainer = std::vector<PartitionTree>;
        
		std::unique_ptr<IntersectionChecker<D, Precision>> intersectionChecker;
        std::variant<ChildContainer, Point_Ids> attachment;
        //size_t numberOfPoints;
        
		PartitionTree() = default;
        
		PartitionTree(const PartitionTree& obj)
			: intersectionChecker(obj.intersectionChecker->copy()), attachment(obj.attachment)
		{
		}

		PartitionTree(PartitionTree&&) = default;

		PartitionTree& operator=(const PartitionTree& rhs)
		{
			intersectionChecker = rhs.intersectionChecker.copy();
			attachment = rhs.attachment;
		}

		PartitionTree& operator=(PartitionTree&&) = default;

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
    std::pair<Partitioning<D, Precision>, std::vector<std::unique_ptr<IntersectionChecker<D, Precision>>>> toOldPartitioning(typename PartitionTree<D, Precision>::ChildContainer children) {
        Partitioning<D, Precision> result;
		std::vector<std::unique_ptr<IntersectionChecker<D, Precision>>> checker;
        size_t id = 0;
        for(auto& childTree : children) {
            childTree.collect();
            auto ids = std::move(std::get<Point_Ids>(childTree.attachment));
            Partition<D, Precision> partition(ids.size());
            partition.points = std::move(ids);
            partition.bounds = childTree.intersectionChecker->bounds();
            partition.id = id++;
            checker.push_back(std::move(childTree.intersectionChecker));
            result.push_back(std::move(partition));
        }
        return std::make_pair(std::move(result), std::move(checker));
    }
}
