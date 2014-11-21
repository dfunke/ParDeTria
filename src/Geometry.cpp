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

dPointStats getPointStats(const dPoints & points, const dPointIds * ids){

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
