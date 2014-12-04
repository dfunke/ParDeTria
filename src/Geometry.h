/*
 * Geometry.h
 *
 *  Created on: Oct 28, 2014
 *      Author: dfunke
 */

#pragma once

#include <iostream>
#include <set>

#include "IndexedVector.hxx"
#include "Logger.h"

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

struct dSphere {
	std::array<tCoordinate,D> center;
	tCoordinate radius;
};

struct dBox {
	std::array<tCoordinate,D> coords;
	tCoordinate dim[D];

	/* tests whether sphere is FULLY contained in box */
	bool contains(const dSphere & sphere) const {

		tCoordinate r2 = sphere.radius * sphere.radius;

		for(uint d = 0; d < D; ++d)
			if(!(coords[d] <= sphere.center[d] && sphere.center[d] <= coords[d] + dim[d]))
				return false;

		//the center of the sphere is within the box
		for(uint d = 0; d < D; ++d){
			auto p = sphere.center;
			//project p to the bounday of box in dimension d closest to center of the sphere
			p[d] = sphere.center[d] < coords[d] + dim[d]/2 ? coords[d] : coords[d] + dim[d];

			tCoordinate dist = 0;
			for(uint i = 0; i < D; ++i)
				dist += (sphere.center[d] - p[d]) * (sphere.center[d] - p[d]);

			if(dist < r2)
				return false;
		}

		return true;

	}
};

//#############################################################################

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

	static constexpr uint cINF = ~(0) ^ ((1 << D) - 1);
};

class dPoints : public IndexedVector<dPoint> {

public:

	dPoints()
		: IndexedVector<dPoint>() { }

	dPoints(const IndexedVector<dPoint> & other)
		: IndexedVector<dPoint>(other) { }

	bool operator==(const std::set<uint> & other) const;

	bool operator!=(const std::set<uint> & other) const {
		return !operator==(other);
	}

};

typedef std::set<uint> Ids;

class Partition {

public:
	uint id;
	Ids points;
	dBox bounds;

};

typedef IndexedVector<Partition> Partitioning;

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

	bool equalVertices(const dSimplex & a) const {
		//compare vertices
		for (uint i = 0; i < D + 1; ++i) {
			bool found = false;
			for (uint j = 0; j < D + 1; ++j) {
				if (vertices[i] == a.vertices[j]) {
					//std::cout << this->vertices[i] << " == " << a.vertices[j];
					found = true;
					break;
				}
			}
			if (!found) {
				//std::cout << false << std::endl;
				return false;
			}
		}

		return true;

	}

	bool equalNeighbors(const dSimplex & a) const {
		return neighbors != a.neighbors;
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

	template <typename Container>
	bool containsAll(const Container & p) const {
		for(const auto & x : p ){
			if(!contains(x))
				return false;
		}

		return true;
	}

	template <typename Container>
	bool containsAny(const Container & p) const {
		for(const auto & x : p ){
			if(contains(x))
				return true;
		}

		return false;
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

	dSphere circumcircle(const dPoints & points) const {

		dSphere sphere;

		tCoordinate D = 2 * (points[vertices[0]].coords[0] * (points[vertices[1]].coords[1] - points[vertices[2]].coords[1])
						   + points[vertices[1]].coords[0] * (points[vertices[2]].coords[1] - points[vertices[0]].coords[1])
						   + points[vertices[2]].coords[0] * (points[vertices[0]].coords[1] - points[vertices[1]].coords[1]));

		sphere.center[0] = (( pow(points[vertices[0]].coords[0],2) + pow(points[vertices[0]].coords[1],2)) * (points[vertices[1]].coords[1] - points[vertices[2]].coords[1])
						   +( pow(points[vertices[1]].coords[0],2) + pow(points[vertices[1]].coords[1],2)) * (points[vertices[2]].coords[1] - points[vertices[0]].coords[1])
						   +( pow(points[vertices[2]].coords[0],2) + pow(points[vertices[2]].coords[1],2)) * (points[vertices[0]].coords[1] - points[vertices[1]].coords[1])) / D;


		sphere.center[1] = (( pow(points[vertices[0]].coords[0],2) + pow(points[vertices[0]].coords[1],2)) * (points[vertices[2]].coords[0] - points[vertices[1]].coords[0])
					  	   +( pow(points[vertices[1]].coords[0],2) + pow(points[vertices[1]].coords[1],2)) * (points[vertices[0]].coords[0] - points[vertices[2]].coords[0])
						   +( pow(points[vertices[2]].coords[0],2) + pow(points[vertices[2]].coords[1],2)) * (points[vertices[1]].coords[0] - points[vertices[0]].coords[0])) / D;

		sphere.radius = sqrt(pow(sphere.center[0] - points[vertices[0]].coords[0], 2) + pow(sphere.center[1] - points[vertices[0]].coords[1], 2));

		return sphere;
	}

	bool isNeighbor(const dSimplex & other) const {
		uint sharedVertices = 0;

		for(uint i = 0; i < D+1; ++i){
			for(uint j = 0; j < D+1; ++j){
				sharedVertices += (vertices[i] == other.vertices[j]);
			}
		}

		assert(sharedVertices <= D || id == other.id);

		//if(sharedVertices > 2 && id != other.id)
		//	PLOG << "Distinct vertices " << *this " and " << other << " share more than 2 hits";

		return sharedVertices == D;
	}

public:
	uint id;
	std::array<uint, D+1> vertices;
	std::set<uint> neighbors;


public:

	static bool isFinite(const uint & i){
		return i != cINF;
	}

	static constexpr uint cINF = ~(0);
};

std::ostream & operator<<(std::ostream & o, const dPoint & p);

std::ostream & operator<<(std::ostream & o, const dSimplex & p);

std::ostream & operator<<(std::ostream & o, const Partition & p);

std::ostream & operator<<(std::ostream & o, const dBox & b);

struct VerificationReport;

class dSimplices : public IndexedVector<dSimplex> {

public:

	dSimplices()
		: IndexedVector<dSimplex>() { }

	dSimplices(const IndexedVector<dSimplex> & other)
		: IndexedVector<dSimplex>(other) { }

	bool operator==(const dSimplices & other) const;

	bool operator!=(const dSimplices & other) const {
		return !operator==(other);
	}

	bool verify(const dPoints & points) const;

	VerificationReport verify(const dSimplices & realDT) const;

	template <typename Container>
	dSimplices findSimplices(const Container & points, const bool all = false) const{

		dSimplices result;

		if (all) {
			for (const auto & s : *this) {
				if (s.containsAll(points))
					result.insert(s);
			}
		} else {
			for (const auto & s : *this) {
				if (s.containsAny(points))
					result.insert(s);
			}
		}

		return result;

	}

};

struct VerificationReport {
	bool valid;
	dSimplices missing;
	dSimplices invalid;
};

struct dPointStats {
	dPoint min;
	dPoint mid;
	dPoint max;
};

dPointStats getPointStats(const dPoints & points, const Ids * = nullptr);

//#############################################################################
