#include "Geometry.h"

#include <algorithm>

std::ostream & operator<<(std::ostream & o, const Point3D & p){
	o << "[" << p.x << ", " << p.y << ", " << p.z << "]";
	return o;
}

std::ostream & operator<<(std::ostream & o, const dPoint & p){
	o << p.id << "-[" << p.coords[0];
	for(uint i = 1; i < D; ++i)
		o << ", " << p.coords[i];
	o << "]";
	return o;
}

std::ostream & operator<<(std::ostream & o, const dSimplex & p){
	o << p.id << ": V = [" << p.vertices[0];
	for(uint i = 1; i < D+1; ++i)
		o << ", " << p.vertices[i];
	o << "] N = [";
	for(auto it = p.neighbors.begin(); it != p.neighbors.end(); ++it){
		if(it != p.neighbors.begin())
			o << ", ";
		o  << *it;
	}
	o << "]";
	return o;
}

dPointStats getPointStats(const dPoints & points, const Ids * ids){

	dPointStats stats;

	dPoints::const_iterator begin, end;
	dPoints projection; //used if needed

	if(ids == nullptr){
		begin = points.begin();
		end = points.end();
	}
	else {
		projection = points.project(*ids);

		begin = projection.begin();
		end = projection.end();
	}

	for(uint dim = 0; dim < D; ++dim){
		auto minmax = std::minmax_element(begin, end,
				[dim] (const dPoint & a, const dPoint & b) {
					return a.coords[dim] < b.coords[dim];
				});

		stats.min.coords[dim] = (*minmax.first).coords[dim];
		stats.max.coords[dim] = (*minmax.second).coords[dim];
		stats.mid.coords[dim] = (stats.max.coords[dim] + stats.min.coords[dim]) / 2;
	}

	return stats;

}

bool dSimplices::operator==(const dSimplices & other) const {

	if(size() != other.size()){
		PLOG << "my size: " << size() << " other size: " << other.size() << std::endl;
		return false;
	}

	for(const auto & otherSimplex : other){

		//find my simplex, compares simplex id or vertices ids
		auto mySimplex = std::find(this->begin(), this->end(), otherSimplex);

		if(mySimplex == this->end()){
			PLOG << "did not find simplex" << otherSimplex << std::endl;
			return false;
		}

		//check neighbors
		if(mySimplex->neighbors != otherSimplex.neighbors){
			//this should never happen
			return false;
		}
	}

	return true;

}

bool dPoints::operator==(const std::set<uint> & other) const {
	if(size() != other.size())
		return false;

	for(const auto & p : other){
		if(!contains(p))
			return false;
	}

	return true;
}

bool dSimplices::verify(const dPoints & points) const{

	INDENT
	bool valid = true;

	//verify that every input point is used
	LOG << "Checking points" << std::endl;
	std::set<uint> usedPoints;
	for(const auto & s : *this){
		usedPoints.insert(s.vertices.begin(), s.vertices.end());
	}
	if(points != usedPoints){
		//not all points of input used
		LOG << "Points of input not used: ";
		for(const auto & p : points)
			if(usedPoints.count(p.id) != 1)
				CONT << p << " ";
		CONT << std::endl;
		LOG << "Used points not in input: ";
		for(const auto & p : usedPoints)
			if(!points.contains(p))
				CONT << p << " ";
		CONT << std::endl;

		valid = false;
	}

	//verify where-used data structure
	LOG << "Checking where-used relation" << std::endl;
	for(const auto & s : *this){
		for(const auto & p : s.vertices){
			//point p of s not correctly flagged as used in s
			if(points[p].simplices.count(s.id) != 1){
				LOG << "Point " << p << " NOT flagged as used in " << s << std::endl;;
				Logger::getInstance().logContainer(points[p].simplices, Logger::Verbosity::NORMAL, "p.simplices");
				valid = false;
			}
		}
	}
	for(const auto & p : points){
		for(const auto & id : p.simplices){

			if(!contains(id))
				continue; //simplex of another triangulation

			//p is flagged as being used in s, but its not
			const auto & s = operator[](id);
			if(std::find(s.vertices.begin(), s.vertices.end(), p.id) == s.vertices.end()){
				LOG << "Point " << p << " SHOULD be used in " << s << std::endl;
				Logger::getInstance().logContainer(p.simplices, Logger::Verbosity::NORMAL, "p.simplices");
				valid = false;
			}
		}
	}

	//verify that all simplices with a shared D-1 simplex are neighbors
	LOG << "Checking neighbors" << std::endl;
	for(const auto & a : *this){
		for(const auto & b : *this){
					//a and b are neighbors: the neighbor property is symmetric and the corresponding simplices must be present in the neighbors arrays accordingly
			if(    ( a.isNeighbor(b) && (!b.isNeighbor(a) || a.neighbors.count(b.id) != 1 || b.neighbors.count(a.id) != 1))
					//a and b are NOT neighbors: must be symmetric and simplices NOT be present in neighbors arrays
				|| (!a.isNeighbor(b) && ( b.isNeighbor(a) || a.neighbors.count(b.id) != 0 || b.neighbors.count(a.id) != 0)) ){

				LOG << "Wrong neighbor relation between " << a << " and " << b << std::endl;
				valid = false;
			}
		}
	}

	//check in circle criterion
	for(const auto & s : *this){
		for(const auto & p : points){
			bool contains = s.contains(p);
			bool inCircle = s.inCircle(p, points);
			if(contains != inCircle){
				LOG << "Point " << p << " is " << (inCircle ? "" : "NOT ") << "in circle of " << s
				    << " but should " << (contains ? "" : "NOT ") << "be" << std::endl;
				valid = false;
			}
		}
	}

	LOG << "Triangulation is " << (valid ? "" : "NOT ") << "valid" << std::endl;

	DEDENT
	return valid;

}
