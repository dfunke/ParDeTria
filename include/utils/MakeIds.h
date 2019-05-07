#pragma once
#include "Geometry.h"

template<uint D, typename Precision>
Point_Ids makePointIds(const dPoints<D, Precision>& points) {
    Point_Ids allPoints(points.size());
    for (tIdType i = 1; i < points.finite_size(); ++i)
        allPoints.insert(points.offset() + i);
    for (tIdType infVertex = 0; infVertex < dPoint<D,Precision>::nINF;
         ++infVertex)

        allPoints.insert(dPoint<D,Precision>::cINF + infVertex);

    return allPoints;
}
