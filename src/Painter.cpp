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

void Painter::draw(const dPoints & points, bool drawPartition){

	for(const auto & p : points){
		draw(p);
	}

	if(drawPartition){
		auto stats = getPointStats(points);

		//draw partition borders
		cr->move_to(translatePoint(stats.mid.coords[0], 0),0);
		cr->line_to(translatePoint(stats.mid.coords[0], 0),height());

		cr->move_to(0, translatePoint(stats.mid.coords[1],1));
		cr->line_to(width(), translatePoint(stats.mid.coords[1],1));

		cr->stroke();
	}
}

void Painter::draw(const dSimplices & simplices, const dPoints & points, bool drawCC){

	for(const auto & s : simplices){
		for(uint d = 0; d < D; ++d){
			line(points, s.vertices[d], s.vertices[d+1]);
		}
		//close the loop
		line(points, s.vertices[D], s.vertices[0]);
	}

	cr->stroke();

	if(drawCC) {
		for(const auto & s : simplices) {
			PLOG << "CC of " << s.id << std::endl;
			INDENT
			circumCircle(s, points);
			DEDENT
		}
	}
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

void Painter::line(const dPoints & points, uint a, uint b){

	auto A = getPoint(points, a);
	auto B = getPoint(points, b);

	if(!A.isFinite() || !B.isFinite()){
		std::vector<double> dash = { 2, 2};
		cr->set_dash(dash, 0);
	}

	cr->move_to(translatePoint(A.coords[0],0), translatePoint(A.coords[1],1));
	cr->line_to(translatePoint(B.coords[0],0), translatePoint(B.coords[1],1));
	cr->stroke();

	cr->unset_dash();

}

void Painter::circumCircle(const dSimplex & s, const dPoints & points){

	dPoint p[D+1];
	for(uint i = 0; i < D+1; ++i){
			p[i] =getPoint(points, s.vertices[i]);
	}

	tCoordinate D = 2 * (p[0].coords[0] * (p[1].coords[1] - p[2].coords[1])
					   + p[1].coords[0] * (p[2].coords[1] - p[0].coords[1])
					   + p[2].coords[0] * (p[0].coords[1] - p[1].coords[1]));

	PLOG << "D = " << D << std::endl;

	dPoint u;
	u.coords[0] = (( pow(p[0].coords[0],2) + pow(p[0].coords[1],2)) * (p[1].coords[1] - p[2].coords[1])
				  +( pow(p[1].coords[0],2) + pow(p[1].coords[1],2)) * (p[2].coords[1] - p[0].coords[1])
				  +( pow(p[2].coords[0],2) + pow(p[2].coords[1],2)) * (p[0].coords[1] - p[1].coords[1])) / D;


	u.coords[1] = (( pow(p[0].coords[0],2) + pow(p[0].coords[1],2)) * (p[2].coords[0] - p[1].coords[0])
				  +( pow(p[1].coords[0],2) + pow(p[1].coords[1],2)) * (p[0].coords[0] - p[2].coords[0])
				  +( pow(p[2].coords[0],2) + pow(p[2].coords[1],2)) * (p[1].coords[0] - p[0].coords[0])) / D;

	//draw circumcenter
	draw(u);

	PLOG <<"u = " << u << std::endl;

	tCoordinate r = sqrt(pow(u.coords[0] - p[0].coords[0], 2) + pow(u.coords[1] - p[0].coords[1], 2));

	PLOG <<"r = " << r << std::endl;

	cr->arc(translatePoint(u.coords[0], 0), translatePoint(u.coords[1],1), translateLength(r,0), 0, 2*M_PI);
	cr->stroke();

}
