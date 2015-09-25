#pragma once

#include "Geometry.h"
#include "Random.h"
#include "SpatialSort.h"

template<uint D, typename Precision>
dPoints<D, Precision> genPoints(const tIdType n, const dBox<D, Precision> &bounds,
                                std::function<Precision()> &dice) {
    dPoints<D, Precision> points;
    points.resize(n+1);

    dVector<D, Precision> dim = bounds.high;
    for (uint d = 0; d < D; ++d) {
        dim[d] -= bounds.low[d];
    }

    for (uint d = 0; d < D; ++d) {
        points[0].coords[d] = 0;
    }

    for (tIdType i = 1; i <= n; ++i) {
        dPoint<D, Precision> p;
        //p.id = i;
        for (uint d = 0; d < D; ++d) {
            p.coords[d] = bounds.low[d] + dim[d] * dice();
        }

        //points.emplace_back(p);
        points[i] = std::move(p);
    }

    CGALSpatialSorter<D, Precision> sort;
    sort.sort(points);

    // TODO checks for colliding points, straight lines, etc.

    return points;
}