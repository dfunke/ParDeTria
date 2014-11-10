
//std library
#include <fstream>
#include <iostream>
#include <random>
#include <chrono>
#include <vector>
#include <functional>
#include <limits>
#include <sstream>
#include <thread>
#include <future>

//root
#include <TGraph2D.h>
#include <TGraphDelaunay2D.h>

//CGAL
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>

//Triangle
#ifndef REAL
	#ifdef SINGLE
		#define REAL float
	#else /* not SINGLE */
		#define REAL double
	#endif /* not SINGLE */
#endif /* ndef REAL */

#include <triangle.h>

//Voronoi
#include <VoronoiDiagramGenerator.h>

//Delaunay
extern "C" {
	#include <delaunay.h>
}

//own
#include "Timings.h"
#include "Random.h"
#include "Geometry.h"

//#############################################################################

typedef std::function<tDuration(const Points3D &)> tTestFunction;

struct tTest {
	tTestFunction fn;
	std::ofstream file;
	std::string name;
};

typedef std::vector<tTest*> tTests;

//#############################################################################

Points3D genPoints(const uint n, const Box3D & bounds, std::function<tCoordinate()> & dice){

	Points3D points(n);


	for(uint i = 0; i < n; ++i){
		points[i].x = bounds.x + bounds.dx * dice();
		points[i].y = bounds.y + bounds.dy * dice();
		points[i].z = bounds.z + bounds.dz * dice();
	}

	//TODO checks for colliding points, straight lines, etc.

	return points;

}

tDurations getTiming(tTestFunction & function, const Points3D & input, const uint reps = 10){

	tDurations meas;
	meas.reserve(reps);


	function(input); //warm-up
	for(uint i = 0; i < reps; ++i)
		meas.push_back(function(input));

	return meas;
}

tDuration delaunayRoot(const Points3D & points){

	TGraph2D graph(points.size());

	for(uint i = 0; i < points.size(); ++i)
		graph.SetPoint(i, points[i].x, points[i].y, points[i].z);

	TGraphDelaunay2D delaunay(&graph);

	auto start = Clock::now();
	delaunay.FindAllTriangles();
	auto end = Clock::now();

	return tDuration(end - start);
}

tDuration delaunayCgal(const Points3D & points){

	typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
	typedef CGAL::Delaunay_triangulation_2<K> Triangulation;

	//copy points into CGAL structure
	std::vector<Triangulation::Point> cPoints;
	cPoints.reserve(points.size());
	for(const auto & p : points)
		cPoints.emplace_back(p.x, p.y);

	Triangulation t;

	auto start = Clock::now();
	t.insert(cPoints.begin(), cPoints.end());
	auto end = Clock::now();

	return tDuration(end - start);


}

tDuration delaunayFortune(const Points3D & points) {


	float * xValues = new float[points.size()];
	float * yValues = new float[points.size()];

	float minX = std::numeric_limits<float>::max();
	float minY = std::numeric_limits<float>::max();

	float maxX = std::numeric_limits<float>::min();
	float maxY = std::numeric_limits<float>::min();

	for(uint i = 0; i < points.size(); ++i){
		xValues[i] = points[i].x;
		yValues[i] = points[i].y;

		minX = std::min(minX, points[i].x);
		maxX = std::max(maxX, points[i].x);

		minY = std::min(minY, points[i].y);
		maxY = std::max(maxY, points[i].y);
	}

	VDG::VoronoiDiagramGenerator vdg;

	auto start = Clock::now();
	vdg.generateVoronoi(xValues,yValues,points.size(), minX, maxX,minY,maxY);
	auto end = Clock::now();

	delete[] xValues;
	delete[] yValues;

	return tDuration(end - start);

}

tDuration delaunayDAC(const Points3D & points) {

	del_point2d_t * dPoints = new del_point2d_t[points.size()];


	for(uint i = 0; i < points.size(); ++i){
		dPoints[i].x = points[i].x;
		dPoints[i].y = points[i].y;
	}

	auto start = Clock::now();
	auto res = delaunay2d_from(dPoints, points.size(), NULL);
	auto end = Clock::now();

	delaunay2d_release(res);
	delete[] dPoints;

	return tDuration(end - start);

}

enum class Algorithms : char {
	FORTUNE = 'F',
	INSERTION = 'i',
	DAC = '\0',
};

tDuration delaunayTriangle(const Points3D & points, Algorithms alg = Algorithms::DAC, bool verbose = false) {

	auto initStruct = [] (triangulateio & s) {
							  s.pointlist = nullptr;                                               /* In / out */
							  s.pointattributelist = nullptr;                                      /* In / out */
							  s.pointmarkerlist = nullptr;                                          /* In / out */
							  s.numberofpoints = 0;                                            /* In / out */
							  s.numberofpointattributes = 0;                                   /* In / out */

							  s.trianglelist = nullptr;                                             /* In / out */
							  s.triangleattributelist = nullptr;                                   /* In / out */
							  s.trianglearealist = nullptr;                                         /* In only */
							  s.neighborlist = nullptr;                                             /* Out only */
							  s.numberoftriangles = 0;                                         /* In / out */
							  s.numberofcorners = 0;                                           /* In / out */
							  s.numberoftriangleattributes = 0;                                /* In / out */

							  s.segmentlist = nullptr;                                              /* In / out */
							  s.segmentmarkerlist = nullptr;                                        /* In / out */
							  s.numberofsegments = 0;                                          /* In / out */

							  s.holelist = nullptr;                        /* In / pointer to array copied out */
							  s.numberofholes = 0;                                      /* In / copied out */

							  s.regionlist = nullptr;                      /* In / pointer to array copied out */
							  s.numberofregions = 0;                                    /* In / copied out */

							  s.edgelist = nullptr;                                                 /* Out only */
							  s.edgemarkerlist = nullptr;            /* Not used with Voronoi diagram; out only */
							  s.normlist = nullptr;                /* Used only with Voronoi diagram; out only */
							  s.numberofedges = 0;                                             /* Out only */
						};

	auto freeStruct = [] (triangulateio & s) {
							  free(s.pointlist);                                               /* In / out */
							  free(s.pointattributelist);                                      /* In / out */
							  free(s.pointmarkerlist);                                          /* In / out */

							  free(s.trianglelist);                                             /* In / out */
							  free(s.triangleattributelist);                                   /* In / out */
							  free(s.trianglearealist);                                         /* In only */
							  free(s.neighborlist);                                             /* Out only */

							  free(s.segmentlist);                                              /* In / out */
							  free(s.segmentmarkerlist);                                        /* In / out */

							  free(s.holelist);                        /* In / pointer to array copied out */

							  free(s.regionlist);                      /* In / pointer to array copied out */

							  free(s.edgelist);                                                 /* Out only */
							  free(s.edgemarkerlist);            /* Not used with Voronoi diagram; out only */
							  free(s.normlist);                /* Used only with Voronoi diagram; out only */
						};

	struct triangulateio in, out, vorout;
	initStruct(in); initStruct(out); initStruct(vorout);

	/* Define input points. */

	in.numberofpoints = points.size();
	in.pointlist = (REAL *) malloc(in.numberofpoints * 2 * sizeof(REAL));

	for (uint i = 0; i < points.size(); ++i) {
		in.pointlist[2 * i] = points[i].x;
		in.pointlist[2 * i + 1] = points[i].y;
	}

	std::stringstream options;
	options << "pczAevn" << ( verbose ? "V" : "Q") << (char) alg;

	auto start = Clock::now();
	triangulate(const_cast<char*>(options.str().c_str()), &in, &out, &vorout);
	auto end = Clock::now();

	freeStruct(in); freeStruct(out); freeStruct(vorout);

	return tDuration(end - start);

}

void runTest(const uint n, const Box3D & bounds, std::function<tCoordinate()> & dice, const uint reps,
		tTests tests){

	//initialize measurement array
	std::vector<tDurations> meas(tests.size());
	for(uint i = 0; i  < tests.size(); ++i){
		meas[i].reserve(reps*reps);
	}

	tDurations tmp;

	for(uint _ = 0; _ < reps; ++_){
		Points3D points = genPoints(n, bounds, dice);

		for(uint i = 0; i < tests.size(); ++i){
			tmp = getTiming(tests[i]->fn, points, reps);
			meas[i].insert(meas[i].end(), tmp.begin(), tmp.end());
		}
	}


	for(uint i = 0; i < tests.size(); ++i){
		std::stringstream ss;

		tMeas stat = stats(meas[i]);
		ss << n  << " " << stat.avg << " " << stat.std << std::endl;

		tests[i]->file << ss.str();
		std::cout << tests[i]->name << " " << ss.str();
	}
}

int main(int argc, char* argv[]) {

	Box3D bounds;
	bounds.x = bounds.y = bounds.z = 0;
	bounds.dx = bounds.dy = bounds.dz = 100;

	std::uniform_real_distribution<tCoordinate> distribution(0,1);
	std::function<tCoordinate()> dice = std::bind (distribution, generator);

	tTests tests;
	std::string header = "n time time_std\n";

	tTest root;
	root.name = "ROOT";
	root.fn = std::bind(delaunayRoot, std::placeholders::_1);
	root.file.open("data/root.csv", std::ofstream::trunc);
	root.file << header;
	tests.push_back(&root);

	tTest cgal;
	cgal.name = "CGAL";
	cgal.fn = std::bind(delaunayCgal, std::placeholders::_1);
	cgal.file.open("data/cgal.csv", std::ofstream::trunc);
	cgal.file << header;
	tests.push_back(&cgal);

	tTest triangleDAC;
	triangleDAC.name = "Triangle_DAC";
	triangleDAC.fn = std::bind(delaunayTriangle, std::placeholders::_1, Algorithms::DAC, false);
	triangleDAC.file.open("data/triangleDAC.csv", std::ofstream::trunc);
	triangleDAC.file << header;
	tests.push_back(&triangleDAC);

	tTest triangleInc;
	triangleInc.name = "Triangle_Inc";
	triangleInc.fn = std::bind(delaunayTriangle, std::placeholders::_1, Algorithms::INSERTION, false);
	triangleInc.file.open("data/triangleInc.csv", std::ofstream::trunc);
	triangleInc.file << header;
	tests.push_back(&triangleInc);

	tTest triangleFortune;
	triangleFortune.name = "Triangle_Fortune";
	triangleFortune.fn = std::bind(delaunayTriangle, std::placeholders::_1, Algorithms::FORTUNE, false);
	triangleFortune.file.open("data/triangleFortune.csv", std::ofstream::trunc);
	triangleFortune.file << header;
	tests.push_back(&triangleFortune);

	tTest fortune;
	fortune.name = "Fortune";
	fortune.fn = std::bind(delaunayFortune, std::placeholders::_1);
	fortune.file.open("data/fortune.csv", std::ofstream::trunc);
	fortune.file << header;
	tests.push_back(&fortune);

	tTest delaunay;
	delaunay.name = "DAC";
	delaunay.fn = std::bind(delaunayDAC, std::placeholders::_1);
	delaunay.file.open("data/dac.csv", std::ofstream::trunc);
	delaunay.file << header;
	tests.push_back(&delaunay);

	std::cout << "Running with " << tests.size() << std::endl;

	uint N = 1e6;
	for(uint n = 10; n < N; n *= 2){

		runTest(n, bounds, dice, 10, tests);

	}


}
