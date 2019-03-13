#pragma once

#include "Partitioner.h"
#include "PartitionTree.h"

namespace LoadBalancing
{ 
    template <uint D, typename Precision, typename MonitorT>
    struct FlatPartitioner : public Partitioner<D, Precision, MonitorT>
    {
	    FlatPartitioner(std::unique_ptr<Partitioner<D, Precision, MonitorT>> originalPartitioner)
            : originalPartitioner(std::move(originalPartitioner))
        {
        }

        PartitionTree<D, Precision> partition(const dBox<D, Precision>& bounds,
                                        dPoints<D, Precision>& points,
                                        const Point_Ids& pointIds,
                                        MonitorT& monitor) override
        {
	        auto tree = originalPartitioner->partition(bounds, points, pointIds, monitor);
	        tree.flatten();
            return tree;
        }
        
        std::string info() const override
        {
            return originalPartitioner->info() + " (flattened)";
        }
        
        const Sampling<D, Precision>& sampling() const {
            return originalPartitioner->sampling();
        }

   private:
		std::unique_ptr<Partitioner<D, Precision, MonitorT>> originalPartitioner;
	};
}

