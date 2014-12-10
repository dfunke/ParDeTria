#pragma once

#include "Geometry.h"

dSimplices delaunayCgal(dPoints &points, const Ids *ids = nullptr,
                        bool filterInfinite = false);
