#pragma once

#include "Geometry.h"

dSimplices delaunayCgal(dPoints & points, const dPointIds * ids = nullptr, bool filterInfinite = false);
