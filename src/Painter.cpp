#include "Painter.h"

template <uint D, typename Precision>
bool Painter<D, Precision>::ENABLED = true;

// specializations

template class Painter<2, float>;
template class Painter<3, float>;

template class Painter<2, double>;
template class Painter<3, double>;
