/*
 * Geometry.h
 *
 *  Created on: Oct 28, 2014
 *      Author: dfunke
 */

#pragma once

#include <vector>
#include <iostream>

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

struct dPoint {
	uint id;
	tCoordinate coords[D];

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

	static inline bool isFinite(const uint & i) {
		//INDENT
		//PLOG << i << " -> " << (i & cINF) << " != " << cINF << ": " << !((i & cINF) == cINF) << std::endl;
		//DEDENT
		return !((i & cINF) == cINF);
	}

	static constexpr uint cINF = ~(0) << 2*D;
};


struct dSimplex {
	uint id;
	uint vertices[D+1];
	uint neighbors[D+1];

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

	bool isFinite() const {
		bool finite = true;

		for(uint i = 0; i < D+1; ++i){
			if(!dPoint::isFinite(vertices[i]))
				finite = false;
		}

		return finite;
	}

	static bool isFinite(const uint & i){
		return i != cINF;
	}

	static constexpr uint cINF = ~(0);
};

std::ostream & operator<<(std::ostream & o, const dPoint & p);

std::ostream & operator<<(std::ostream & o, const dSimplex & p);

typedef std::vector<dPoint> dPoints;
typedef std::vector<dPoints> Partition;
typedef std::vector<dSimplex> dSimplices;

struct dBox {
	tCoordinate coords[D];
	tCoordinate dim[D];
};

struct dPointStats {
	dPoint min;
	dPoint mid;
	dPoint max;
};

dPointStats getPointStats(const dPoints & points);

//#############################################################################
