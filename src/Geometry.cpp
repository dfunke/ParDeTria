#include "Geometry.h"

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
	o << "] N = [" << p.neighbors[0];
	for(uint i = 1; i < D+1; ++i)
		o << ", " << p.neighbors[i];
	o << "]";
	return o;
}

dPointStats getPointStats(const dPoints & points, const dPointIds * ids){

	dPointStats stats;

	dPoints::const_iterator begin, end;
	dPoints projection;

	if(ids == nullptr){
		begin = points.begin();
		end = points.end();
	}
	else {
		projection.reserve(ids->size());
		for(const auto & i : *ids)
			projection.push_back(points[i]);

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
