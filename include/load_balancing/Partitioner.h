#pragma once
#include <string>
#include "Geometry.h"
#include "PartitionTree.h"

namespace LoadBalancing
{
    template <uint D, typename Precision>
    struct Partitioner
    {
        virtual PartitionTree<D, Precision> partition(const dBox<D, Precision>& bounds,
                                        const dPoints<D, Precision>& points,
                                        const Point_Ids& pointIds) = 0;
        
        virtual std::string info() const = 0;
    };
}
