#include <functional>
#include <fstream>
#include <iostream>

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

double sinsin(double x, double y){
	return (sin(x) / x) * (sin(y) / y) + 0.2;
}

//******************************************************************************

Points3D genSample(const uint n, const Box & bounds,
		tTestFunction & f,
		std::function<tCoordinate()> & dice){

	Points3D points(n);

	tCoordinate pointsPerDim = std::sqrt(n);
	tCoordinate stepX = bounds.dx / pointsPerDim;
	tCoordinate stepY = bounds.dy / pointsPerDim;

	for(uint i = 0; i < pointsPerDim; ++i){
		for(uint j = 0; j < pointsPerDim; ++j){
			points[i].x = bounds.x + i * stepX + stepX*(dice() - .5);
			if(points[i].x < bounds.x) points[i].x = bounds.x;
			if(points[i].x > bounds.x + bounds.dx) points[i].x =  bounds.x + bounds.dx;

			points[i].y = bounds.y + j * stepY + stepY*(dice() - .5);
			if(points[i].y < bounds.y) points[i].y = bounds.y;
			if(points[i].y > bounds.y + bounds.dy) points[i].y =  bounds.y + bounds.dy;

			points[i].z = f(points[i].x, points[i].y);
		}
	}

	return points;

}

//******************************************************************************

int main(int argc, char* argv[]) {

	std::ofstream timings("timings_" SUFFIX ".csv", std::ofstream::trunc);
	timings << "#n gen gen_std ip ip_std" << std::endl;

	std::ofstream precision("precision_" SUFFIX ".csv", std::ofstream::trunc);
	precision << "#n x y z ip" << std::endl;

	uint N = 1e5; //maximum number of measurement points
	uint R = 10; //number of repetitions
	uint IP = 1e4; //number of interpolation points

	Box bounds;
	bounds.x = bounds.y = -6; bounds.z = 0;
	bounds.dx = bounds.dy = 12; bounds.dz = 0;

	std::uniform_real_distribution<tCoordinate> distribution(0,1);
	std::function<tCoordinate()> dice = std::bind (distribution, generator);

	tTestFunction testFunction = std::bind(sinsin, std::placeholders::_1, std::placeholders::_2);

	std::cout << "Generating " << N << " sample points" << std::endl;
	Points3D samplePoints = genSample(N, bounds, testFunction, dice);

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

			auto start = Clock::now();
			tCoordinate ip = dt.Eval(samplePoints[i].x, samplePoints[i].y);
			auto end = Clock::now();

			meas.emplace_back(end - start);

			precision << n << " " << samplePoints[i].x << " " << samplePoints[i].y << " " << samplePoints[i].z << " " << ip << std::endl;
		}

		auto ipStats = stats(meas);
		timings << n << " " << genStats.avg << " " << genStats.std
				<< " " << ipStats.avg << " " << ipStats.std << std::endl;
	}

	timings.close();
	precision.close();
}
