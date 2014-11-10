#include <functional>
#include <fstream>
#include <iostream>
#include <algorithm>

//root
#include <TGraph2D.h>
#include <TGraphDelaunay2D.h>

//own
#include "Timings.h"
#include "Random.h"
#include "Geometry.h"

#ifdef HAS_CGAL
	#define SUFFIX "cgal"
#else
	#define SUFFIX "tria"
#endif

//******************************************************************************
// test functions

typedef std::function<double(double, double)> tTestFunction;
typedef std::vector<tTestFunction> tTestFunctions;

double fSin(double x, double y){
	return sin(4*x*y);
}

double fSinSin(double x, double y){
	return (sin(x) / x) * (sin(y) / y) + 0.2;
}

double fCubic(double x, double y){
	return pow(x,3)-3*x+pow(y,3)-3*y;
}

double fSqrtSinCos(double x, double y){
	return sqrt(sin(x)+cos(y));
}

//******************************************************************************

Points3D genSample(const uint n, const Box3D & bounds,
		tTestFunction & f,
		std::function<tCoordinate()> & dice){

	Points3D points(n);
	uint idx = 0;

	tCoordinate pointsPerDim = std::sqrt(n);
	tCoordinate stepX = bounds.dx / pointsPerDim;
	tCoordinate stepY = bounds.dy / pointsPerDim;

	for(uint i = 0; i < std::floor(pointsPerDim); ++i){
		for(uint j = 0; j < std::floor(pointsPerDim); ++j){
			points[idx].x = bounds.x + i * stepX + stepX*(dice() - .5);
			if(points[idx].x < bounds.x) points[idx].x = bounds.x;
			if(points[idx].x > bounds.x + bounds.dx) points[idx].x =  bounds.x + bounds.dx;

			points[idx].y = bounds.y + j * stepY + stepY*(dice() - .5);
			if(points[idx].y < bounds.y) points[idx].y = bounds.y;
			if(points[idx].y > bounds.y + bounds.dy) points[idx].y =  bounds.y + bounds.dy;

			points[idx].z = f(points[idx].x, points[idx].y);

//			std::cout << idx << ":" << points[idx] << std::endl;

			++idx;
		}
	}

	//this might lead to less points than n
	while(idx < n){
		points[idx].x = bounds.x + bounds.dx*dice();
		points[idx].y = bounds.y + bounds.dy*dice();
		points[idx].z = f(points[idx].x, points[idx].y);

		++idx;
	}

	//sort the points again lexicographically
	std::sort(points.begin(), points.end(),
		[] (const Point3D & a, const Point3D & b){
			return a.x < b.x && a.y < b.y;
		});

	return points;

}

//******************************************************************************

int main(int argc, char* argv[]) {

	std::ofstream timings("timings_" SUFFIX ".csv", std::ofstream::trunc);
	timings << "#f n gen gen_std ip ip_std" << std::endl;

	std::ofstream precision("precision_" SUFFIX ".csv", std::ofstream::trunc);
	precision << "#f n x y z ip" << std::endl;

	uint N = 1e5; //maximum number of measurement points
	uint R = 10; //number of repetitions
	uint IP = 1e4; //number of interpolation points

	Box3D bounds;
	bounds.x = bounds.y = -6; bounds.z = 0;
	bounds.dx = bounds.dy = 12; bounds.dz = 0;

	std::uniform_real_distribution<tCoordinate> distribution(0,1);
	std::function<tCoordinate()> dice = std::bind (distribution, generator);

	tTestFunctions testFunctions;
	testFunctions.push_back(std::bind(fSin, std::placeholders::_1, std::placeholders::_2));
	testFunctions.push_back(std::bind(fSinSin, std::placeholders::_1, std::placeholders::_2));
	testFunctions.push_back(std::bind(fCubic, std::placeholders::_1, std::placeholders::_2));
	testFunctions.push_back(std::bind(fSqrtSinCos, std::placeholders::_1, std::placeholders::_2));

	for(uint f = 0; f < testFunctions.size(); ++f){

		std::cout << "Generating " << N << " sample points for function " << f << std::endl;
		Points3D samplePoints = genSample(N, bounds, testFunctions[f], dice);

		for (uint n = 10; n < N; n *= 2) {

			std::cout << "Testing with " << n << " points" << std::endl;

			std::cout << "\tFilling graph" << std::endl;
			TGraph2D graph(n);

			double sampleStep = samplePoints.size() / n;
			std::uniform_int_distribution<int> distributionT(-sampleStep/2, sampleStep/2);
			auto diceT = std::bind(distributionT, generator);

			for (uint i = 0; i < n; ++i) {

				int idx = i * sampleStep + diceT();
				if(idx < 0) idx = 0;
				if(idx >= samplePoints.size() - 1) idx = samplePoints.size() - 1;

	//			std::cout << i << "->" << idx << ": " << samplePoints[idx] << std::endl;

				graph.SetPoint(i, samplePoints[idx].x, samplePoints[idx].y, samplePoints[idx].z);
			}

			std::cout << "\tRunning timing tests" << std::endl;
			tDurations meas;
			for (uint r = 0; r < R; ++r) {
				TGraphDelaunay2D dt(&graph);

				auto start = Clock::now();
				dt.FindAllTriangles();
				auto end = Clock::now();

				meas.emplace_back(end - start);
			}

			auto genStats = stats(meas);

			std::cout << "\tRunning precision tests" << std::endl;
			TGraphDelaunay2D dt(&graph);
			dt.FindAllTriangles();

			meas.clear();

			sampleStep = samplePoints.size() / IP;
			std::uniform_int_distribution<int> distributionIP(-sampleStep/2, sampleStep/2);
			auto diceIP = std::bind(distributionIP, generator);

			for(uint i = 0; i < IP; ++i) {

				int idx = i * sampleStep + diceIP();
				if(idx < 0) idx = 0;
				if(idx >= samplePoints.size() - 1) idx = samplePoints.size() - 1;

	//			std::cout << i << "->" << idx << ": " << samplePoints[idx] << std::endl;

				auto start = Clock::now();
				tCoordinate ip = dt.Eval(samplePoints[idx].x, samplePoints[idx].y);
				auto end = Clock::now();

				meas.emplace_back(end - start);

				precision << f << " " << n << " "
						  << samplePoints[idx].x << " " << samplePoints[idx].y << " " << samplePoints[idx].z << " "
						  << ip << std::endl;
			}

			auto ipStats = stats(meas);
			timings << f << " " << n << " " << genStats.avg << " " << genStats.std
					<< " " << ipStats.avg << " " << ipStats.std << std::endl;
		}
	}

	timings.close();
	precision.close();
}
