/*
 * Random.h
 *
 *  Created on: Oct 28, 2014
 *      Author: dfunke
 */

#pragma once

#include <random>

#include "Geometry.h"

//#############################################################################

const uint SEED = 1986;
std::random_device rd;
std::mt19937 generator(SEED);

/* Prototype
 * std::uniform_int_distribution<uint> distribution(0,100);
 * auto dice = std::bind ( distribution, generator );
 */

dPoints genPoints(const uint n, const dBox &bounds,
                  std::function<tCoordinate()> &dice) {
  dPoints points;

  for (uint i = 0; i < n; ++i) {
    points[i].id = i;
    for (uint d = 0; d < D; ++d) {
      points[i].coords[d] = bounds.coords[d] + bounds.dim[d] * dice();
    }
  }

  // TODO checks for colliding points, straight lines, etc.

  return points;
};
