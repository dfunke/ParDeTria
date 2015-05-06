//
// Created by dfunke on 4/8/15.
//

#include "Triangulator.h"

template<uint D, typename Precision>
constexpr char Triangulator<D, Precision>::TOP;

template<uint D, typename Precision>
Ids Triangulator<D, Precision>::allPoints() const {
    Ids allPoints;
    allPoints.reserve(points.size());
    for (uint i = 0; i < points.finite_size(); ++i)
        allPoints.insert(i);
    for (uint infVertex = points.finite_size(); infVertex < points.size();
         ++infVertex)

        allPoints.insert(points[infVertex].id);

    return allPoints;
}

// specializations
template
class Triangulator<2, float>;

template
class Triangulator<3, float>;

template
class Triangulator<2, double>;

template
class Triangulator<3, double>;