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

template <typename Precision>
class DistributionFactory{

public:
  static std::function<Precision()> make(const unsigned char type){

    switch(type){
      case 'u':
      default:
        std::uniform_real_distribution<Precision> distribution(0, 1);
        return std::bind(distribution, generator);
    }
  }

};

template <uint D, typename Precision>
dPoints<D, Precision> genPoints(const uint n, const dBox<D, Precision> &bounds,
                                std::function<Precision()> &dice) {
  dPoints<D, Precision> points;
  points.reserve(n);

  dVector<D, Precision> dim = bounds.high;
  for (uint d = 0; d < D; ++d) {
    dim[d] -= bounds.low[d];
  }

  for (uint i = 0; i < n; ++i) {
    dPoint<D, Precision> p;
    p.id = i;
    for (uint d = 0; d < D; ++d) {
      p.coords[d] = bounds.low[d] + dim[d] * dice();
    }

    points.emplace_back(p);
  }

  // TODO checks for colliding points, straight lines, etc.

  return points;
}
