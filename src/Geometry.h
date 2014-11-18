/*
 * Geometry.h
 *
 *  Created on: Oct 28, 2014
 *      Author: dfunke
 */

#pragma once

#include <iostream>
#include <set>

#include "Logger.h"
#include "IdVector.hxx"

typedef float tCoordinate;

struct Point3D {
	tCoordinate x,y,z;
};

typedef std::vector<Point3D> Points3D;

struct Box3D {
	tCoordinate x,y,z;
	tCoordinate dx,dy,dz;
};

std::ostream & operator<<(std::ostream & o, const Point3D & p);

//#############################################################################
// forward declarations
//
//struct Point;
//struct dSimplex;
//
//std::ostream & operator<<(std::ostream & o, const Point & p);
//std::ostream & operator<<(std::ostream & o, const dSimplex & p);
//
//#############################################################################

//dimensionality of our problem
const uint D = 2;

class dPoint {

public:

	bool operator==(const dPoint & a) const {

		//COUT << "Comparing POINTS THIS " << *this << " and OTHER " << a << ": ";

		if(!(isFinite() ^ a.isFinite())){
			//either none or both points are infinity
			//compare ids
			//std::cout << (this->id == a.id) << std::endl;
			return id == a.id;
		}
		else{
			//either one is infinity
			//compare coordinates
			for(uint i = 0; i < D; ++i){
				if(coords[i] != a.coords[i]){
					//std::cout << false << std::endl;
					return false;
				}
			}
			//std::cout << true << std::endl;
			return true;
		}
	}

	bool operator==(const uint & a) const {
		return id == a;
	}

	inline bool isFinite() const {
		return isFinite(id);
	}

public:

	uint id;
	std::array<tCoordinate, D> coords;
	std::set<uint> simplices;

public:

	static inline bool isFinite(const uint & i) {
		//INDENT
		//PLOG << i << " -> " << (i & cINF) << " != " << cINF << ": " << !((i & cINF) == cINF) << std::endl;
		//DEDENT
		return !((i & cINF) == cINF);
	}

	static constexpr uint cINF = ~(0) << 2*D;
};

typedef IdVector<dPoint> dPoints;
typedef std::vector<uint> dPointIds;
typedef std::vector<dPointIds> Partition;

class dSimplex {

public:

	bool operator==(const dSimplex & a) const {

		//COUT << "Comparing SIMPLICES THIS " << *this << " and OTHER " << a << ": ";

		if(!((id == cINF) ^ (a.id == cINF))){
			//either none or both simplices are infinity
			//compare ids
			//std::cout << (this->id == a.id) << std::endl;
			return id == a.id;
		}
		else {
			//either one is infinity
			//compare vertices
			for(uint i = 0; i < D+1; ++i){
				bool found = false;
				for(uint j = 0; j < D+1; ++j){
					if(vertices[i] == a.vertices[j]){
						//std::cout << this->vertices[i] << " == " << a.vertices[j];
						found = true;
						break;
					}
				}
				if(!found){
					//std::cout << false << std::endl;
					return false;
				}
			}
			//std::cout << true << std::endl;
			return true;
		}
	}

	bool operator==(const uint & a) const {
		return id == a;
	}

	bool contains(uint p) const {
		for (uint d = 0; d < D + 1; ++d) {
			if (p == vertices[d])
				return true;
		}

		return false;
	}

	bool contains(const dPoint & p) const {
		for(uint d = 0; d < D+1; ++d){
			if(p == vertices[d])
				return true;
		}

		return false;
	}

	bool contains(const dPoints & points) const {
		for(const auto & p : points){
			if(!contains(p))
				return false;
		}

		return true;
	}

	bool isFinite() const {
		bool finite = true;

		for(uint i = 0; i < D+1; ++i){
			if(!dPoint::isFinite(vertices[i]))
				finite = false;
		}

		return finite;
	}

	/*
	 * return > 0: abc are counter-clockwise
	 * return < 0: abc clockwise
	 */
	tCoordinate orientation(const dPoints & points) const {
	  tCoordinate acx, bcx, acy, bcy;

	  acx = points[vertices[0]].coords[0] - points[vertices[2]].coords[0];
	  bcx = points[vertices[1]].coords[0] - points[vertices[2]].coords[0];
	  acy = points[vertices[0]].coords[1] - points[vertices[2]].coords[1];
	  bcy = points[vertices[1]].coords[1] - points[vertices[2]].coords[1];
	  return acx * bcy - acy * bcx;
	}

	bool inCircle(const dPoint & p, const dPoints & points) const {
	  tCoordinate adx, ady, bdx, bdy, cdx, cdy;
	  tCoordinate abdet, bcdet, cadet;
	  tCoordinate alift, blift, clift;
	  tCoordinate det;

	  adx = points[vertices[0]].coords[0] - p.coords[0];
	  ady = points[vertices[0]].coords[1] - p.coords[1];
	  bdx = points[vertices[1]].coords[0] - p.coords[0];
	  bdy = points[vertices[1]].coords[1] - p.coords[1];
	  cdx = points[vertices[2]].coords[0] - p.coords[0];
	  cdy = points[vertices[2]].coords[1] - p.coords[1];

	  abdet = adx * bdy - bdx * ady;
	  bcdet = bdx * cdy - cdx * bdy;
	  cadet = cdx * ady - adx * cdy;
	  alift = adx * adx + ady * ady;
	  blift = bdx * bdx + bdy * bdy;
	  clift = cdx * cdx + cdy * cdy;

	  /*
	   * det > 0: d inside  abc (abc counter-clockwise)
	   * 		  d outside abc (abc clockwise)
	   * det < 0: d outside abc (abc counter-clockwise)
	   * 		  d inside  abc (abc clockwise)
	   */

	  det = alift * bcdet + blift * cadet + clift * abdet;

	  return orientation(points) * det >= 0;
	}

public:
	uint id;
	std::array<uint, D+1> vertices;
	std::array<uint, D+1> neighbors;


public:

	static bool isFinite(const uint & i){
		return i != cINF;
	}

	static constexpr uint cINF = ~(0);
};

std::ostream & operator<<(std::ostream & o, const dPoint & p);

std::ostream & operator<<(std::ostream & o, const dSimplex & p);

typedef IdVector<dSimplex> dSimplices;

struct dBox {
	tCoordinate coords[D];
	tCoordinate dim[D];
};

struct dPointStats {
	dPoint min;
	dPoint mid;
	dPoint max;
};

dPointStats getPointStats(const dPoints & points, const dPointIds * = nullptr);

//#############################################################################
