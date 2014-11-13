#include "Painter.h"


void Painter::draw(const dPoint & point){

	if(!point.isFinite()){
		cr->save();
		cr->set_source_rgb(1,1,0);
	}

	cr->arc(translatePoint(point.coords[0], 0), translatePoint(point.coords[1],1), 5, 0, 2*M_PI);
	cr->fill();

	if(!point.isFinite())
		cr->restore();
}

void Painter::draw(const dPoints & points){

	for(const auto & p : points)
		draw(p);
}

void Painter::draw(const dSimplex & simplex, const dPoints & points, bool drawInfinite){

	for(uint d = 0; d < D; ++d){
		line(points, simplex.vertices[d], simplex.vertices[d+1], drawInfinite);
	}
	//close the loop
	line(points, simplex.vertices[D], simplex.vertices[0], drawInfinite);

	cr->stroke();
}


void Painter::draw(const dSimplices & simplices, const dPoints & points, bool drawInfinite){

	for(const auto & s : simplices)
		draw(s, points, drawInfinite);
}

void Painter::drawCircumCircle(const dSimplex & s, const dPoints & points, bool drawInfinite){

	dPoint p[D+1];
	bool allFinite = true;

	for(uint i = 0; i < D+1; ++i){
			p[i] =getPoint(points, s.vertices[i]);
			if(!p[i].isFinite())
				allFinite = false;
	}

	if(!drawInfinite && !allFinite)
		return;

	tCoordinate D = 2 * (p[0].coords[0] * (p[1].coords[1] - p[2].coords[1])
					   + p[1].coords[0] * (p[2].coords[1] - p[0].coords[1])
					   + p[2].coords[0] * (p[0].coords[1] - p[1].coords[1]));

	dPoint u;
	u.coords[0] = (( pow(p[0].coords[0],2) + pow(p[0].coords[1],2)) * (p[1].coords[1] - p[2].coords[1])
				  +( pow(p[1].coords[0],2) + pow(p[1].coords[1],2)) * (p[2].coords[1] - p[0].coords[1])
				  +( pow(p[2].coords[0],2) + pow(p[2].coords[1],2)) * (p[0].coords[1] - p[1].coords[1])) / D;


	u.coords[1] = (( pow(p[0].coords[0],2) + pow(p[0].coords[1],2)) * (p[2].coords[0] - p[1].coords[0])
				  +( pow(p[1].coords[0],2) + pow(p[1].coords[1],2)) * (p[0].coords[0] - p[2].coords[0])
				  +( pow(p[2].coords[0],2) + pow(p[2].coords[1],2)) * (p[1].coords[0] - p[0].coords[0])) / D;

	//draw circumcenter
	draw(u);

	tCoordinate r = sqrt(pow(u.coords[0] - p[0].coords[0], 2) + pow(u.coords[1] - p[0].coords[1], 2));

	PLOG <<"r = " << r << std::endl;

	cr->arc(translatePoint(u.coords[0], 0), translatePoint(u.coords[1],1), translateLength(r,0), 0, 2*M_PI);
	cr->stroke();

}

void Painter::drawCircumCircle(const dSimplices & s, const dPoints & points, bool drawInfinite){
	for(const auto & x : s)
		drawCircumCircle(x,points, drawInfinite);
}

void Painter::drawPartition(const dPoints & points){
	auto stats = getPointStats(points);

	//draw partition borders
	cr->move_to(translatePoint(stats.mid.coords[0], 0),0);
	cr->line_to(translatePoint(stats.mid.coords[0], 0),height());

	cr->move_to(0, translatePoint(stats.mid.coords[1],1));
	cr->line_to(width(), translatePoint(stats.mid.coords[1],1));

	cr->stroke();
}

void Painter::drawNeighbors(const dSimplex & simplex, const dSimplices & neighbors, const dPoints & points, bool drawInfinite){
	for(uint n = 0; n < D+1; ++n){
		if(dSimplex::isFinite(simplex.neighbors[n])){
			draw(neighbors[simplex.neighbors[n]], points, drawInfinite);
		}
	}
}

void Painter::drawNeighbors(const dSimplices & simplices, const dSimplices & neighbors, const dPoints & points, bool drawInfinite){
	for(const auto & x : simplices)
		drawNeighbors(x,neighbors, points, drawInfinite);
}

void Painter::_init(const dBox & _bounds, uint _resolution) {
	bounds = _bounds;

	for(uint d = 0; d < D; ++d){
		offset.coords[d] = bounds.dim[d]; //offset bounds in middle of image

		img.coords[d] = 0; //image starts at 0
		img.dim[d] = 3*bounds.dim[d] * _resolution; //9 quadrants
	}

	cs = Cairo::ImageSurface::create(Cairo::FORMAT_ARGB32,
			width(), height());
	cr = Cairo::Context::create(cs);

	//draw background white
	cr->set_source_rgb(1,1,1);
	cr->paint();

	cr->set_line_width(1.0);
	cr->set_source_rgb(0,0,0);

	//draw bounds
	cr->rectangle(translatePoint(bounds.coords[0], 0),
				  translatePoint(bounds.coords[1], 1),
				  translateLength(bounds.dim[0], 0),
				  translateLength(bounds.dim[1], 1));
	cr->stroke();
}

void Painter::_copy(const Painter & a) {
	bounds = a.bounds;
	img = a.img;
	offset = a.offset;

	cs = Cairo::ImageSurface::create(Cairo::FORMAT_ARGB32,
			width(), height());
	cr = Cairo::Context::create(cs);

	//copy a's surface
	cr->save();
	cr->set_source(a.cs, 0, 0);
	cr->paint();
	cr->restore();
}

const dPoint & Painter::getPoint(const dPoints & points, uint a) const {

	if(dPoint::isFinite(a)){
		return points[a];
	} else {
		auto p = std::find(points.rbegin(), points.rend(), a);
		if(p != points.rend())
			return *p;
		else
			VLOG << "Could not find point " << a << std::endl;
	}

	return dPoint();

}

void Painter::line(const dPoints & points, uint a, uint b, bool drawInfinite){

	auto A = getPoint(points, a);
	auto B = getPoint(points, b);

	if(!A.isFinite() || !B.isFinite()){
		if(!drawInfinite)
			return;

		std::vector<double> dash = { 2, 2};
		cr->set_dash(dash, 0);
	}

	cr->move_to(translatePoint(A.coords[0],0), translatePoint(A.coords[1],1));
	cr->line_to(translatePoint(B.coords[0],0), translatePoint(B.coords[1],1));
	cr->stroke();

	cr->unset_dash();

}
