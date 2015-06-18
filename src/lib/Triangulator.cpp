//
// Created by dfunke on 4/8/15.
//

#include "Triangulator.h"

template<uint D, typename Precision>
constexpr char Triangulator<D, Precision>::TOP;

template<uint D, typename Precision>
Ids Triangulator<D, Precision>::allPoints() const {
    Ids allPoints(points.size());
    for (uint i = 1; i < points.finite_size(); ++i)
        allPoints.insert(i);
    for (uint infVertex = 0; infVertex < dPoint<D,Precision>::nINF;
         ++infVertex)

        allPoints.insert(dPoint<D,Precision>::cINF + infVertex);

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