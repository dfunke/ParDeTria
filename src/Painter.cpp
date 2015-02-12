#include "Painter.h"

#include "Painter_2.cxx"
#include "Painter_3.cxx"

template <uint D, typename Precision>
bool Painter<D, Precision>::ENABLED = true;

// specializations

template class Painter<2, float>;
template class Painter<3, float>;
