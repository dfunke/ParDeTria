#pragma once

//stl
#include <functional>

// own
#include "Geometry.h"

class Triangulator{

public:

	Triangulator(dBox & _bounds, dPoints & _points)
		: bounds(_bounds), points(_points) { }

	Triangulator(const uint n, dBox & _bounds,
				 std::function<tCoordinate()> &dice)
		: bounds(_bounds) {

		points = genPoints(n, dice);
	}

	dSimplices triangulate();

private:
	dSimplices triangulateBase(const Ids partitionPoints,
			const std::string provenance);

	dSimplices triangulateDAC(const Ids partitionPoints,
			const std::string provenance);

	dPoints genPoints(const uint n, std::function<tCoordinate()> &dice);

	Partitioning partition(const Ids &ids);

	Ids getEdge(const dSimplices &simplices, const Partition &partition);

	Ids extractPoints(const Ids &edgeSimplices, const dSimplices &simplices);

	void updateNeighbors(dSimplices &simplices, const std::string &provenance);

	void eliminateDuplicates(dSimplices &DT);

	dSimplices mergeTriangulation(std::vector<dSimplices> &partialDTs,
			const Ids &edgeSimplices, const dSimplices &edgeDT,
			const Partitioning &partitioning, const std::string &provenance,
			const dSimplices *realDT = nullptr);

	void evaluateVerificationReport(const VerificationReport &vr,
			const std::string &provenance) const;

	void evaluateCrossCheckReport(const CrossCheckReport &ccr,
			const std::string &provenance, const dSimplices &DT,
			const dSimplices &realDT) const;

private:
	dBox bounds;
	dPoints points;

	static constexpr tCoordinate SAFETY = 100;
	static const uint BASE_CASE = 100;

};
