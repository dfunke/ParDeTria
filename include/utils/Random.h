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

template<uint D, typename Precision>
class RandomPointHolder {

public:
  typedef std::pair<unsigned char, uint> point_key;
  typedef std::map<point_key, dPoints<D, Precision>> point_map;

  void fill(unsigned char dist, uint N, const dBox<D, Precision> &bounds){

    auto dice = DistributionFactory<Precision>::make(dist);

    //loop over number of points
    uint nPointsIterations = 9 * (log10(N) - 1) + 1;
    for (uint i = 0; i < nPointsIterations; ++i) {
      uint nPoints = ((i % 9) + 1) * std::pow(10, std::floor(i / 9 + 1));

      auto points = genPoints<D, Precision>(nPoints, bounds, dice);

      m_points.emplace(point_key(dist, nPoints), points);

    }
  }

  dPoints<D, Precision> get(const point_key & key) const {
    return dPoints<D,Precision>(m_points.at(key));
  }

  dPoints<D, Precision> get(const unsigned char dist, const uint n) const {
    return get(point_key(dist, n));
  }

private:
  point_map m_points;

};
