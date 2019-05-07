//
// Created by dfunke on 4/8/15.
//

#include "Triangulator.h"
#include "utils/MakeIds.h"

template<uint D, typename Precision>
constexpr char Triangulator<D, Precision>::TOP;

template<uint D, typename Precision>
Point_Ids Triangulator<D, Precision>::allPoints() const {
    return  makePointIds(points);
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
