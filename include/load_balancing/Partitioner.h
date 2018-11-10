#pragma once
#include <string>
#include "Geometry.h"
#include "PartitionTree.h"

namespace LoadBalancing
{
    template <uint D, typename Precision, typename MonitorT>
    struct Partitioner
    {
        virtual PartitionTree<D, Precision> partition(const dBox<D, Precision>& bounds,
                                        dPoints<D, Precision>& points,
                                        const Point_Ids& pointIds,
                                        MonitorT& monitor) = 0;
        
        virtual std::string info() const = 0;
    };
}
