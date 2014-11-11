#pragma once

#include <cairo/cairo.h>
#include <cairomm/cairomm.h>

#include "Geometry.h"
#include "Logger.h"

class Painter {

public:

	Painter(const dBox & _bounds, const dBox & _img)
	: bounds(_bounds), img(_img) {

		cs = Cairo::ImageSurface::create(Cairo::FORMAT_ARGB32,
				width(), height());
		cr = Cairo::Context::create(cs);

		//draw background white
		cr->set_source_rgb(1,1,1);
		cr->paint();

		cr->set_line_width(1.0);
		cr->set_source_rgb(0,0,0);
	}

	Painter(const Painter & a)
	: bounds(a.bounds), img(a.img) {

		cs = Cairo::ImageSurface::create(Cairo::FORMAT_ARGB32,
				width(), height());
		cr = Cairo::Context::create(cs);

		//copy a's surface
		cr->set_source(a.cs, 0, 0);
		cr->paint ();
	}

	void draw(const dPoint & point){
		cr->arc(translate(point.coords[0], 0), translate(point.coords[1],1), 5, 0, 2*M_PI);
		cr->fill();
	}

	void draw(const dPoints & points, bool drawPartition = true){

		for(const auto & p : points){
			draw(p);
		}

		if(drawPartition){
			dPoint midpoint;
			for(uint dim = 0; dim < D; ++dim){
				auto minmax = std::minmax_element(points.begin(), points.end(),
						[dim] (const dPoint & a, const dPoint & b) {
							return a.coords[dim] < b.coords[dim];
						});
				midpoint.coords[dim] = ((*minmax.second).coords[dim] - (*minmax.first).coords[dim]) / 2;
			}

			//draw partition borders
			cr->move_to(translate(midpoint.coords[0], 0),0);
			cr->line_to(translate(midpoint.coords[0], 0),height());

			cr->move_to(0, translate(midpoint.coords[1],1));
			cr->line_to(width(), translate(midpoint.coords[1],1));

			cr->stroke();
		}
	}

	void draw(const dSimplices & simplices, const dPoints & points, bool drawCC = false){

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

	void setColor(tCoordinate r, tCoordinate g, tCoordinate b){
		if(r > 1 || g > 1 || b > 1){
			//rescale colors
			r /= 255; g /= 255; b /= 255;
		}

		cr->set_source_rgb(r, g, b);
	}

	void setLineWidth(tCoordinate w){
		cr->set_line_width(w);
	}

	void savePNG(const std::string & file) const {
		cs->write_to_png(file);
	}

private:

	void line(const dPoints & points, uint a, uint b){
		if(!(dPoint::isFinite(a) && dPoint::isFinite(b)))
			return;

		cr->move_to(translate(points[a].coords[0],0), translate(points[a].coords[1],1));
		cr->line_to(translate(points[b].coords[0],0), translate(points[b].coords[1],1));
	}

	void circumCircle(const dSimplex & s, const dPoints & points){

		dPoint p[D+1];
		for(uint i = 0; i < D+1; ++i){
			if(!dPoint::isFinite(s.vertices[i]))
				return;
			else
				p[i] = points[s.vertices[i]];
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

		cr->arc(translate(u.coords[0], 0), translate(u.coords[1],1), translate(r,0), 0, 2*M_PI);
		cr->stroke();

	}

	tCoordinate width() const {
		return img.dim[0] - img.coords[0];
	}

	tCoordinate height() const {
		return img.dim[1] - img.coords[1];
	}

	tCoordinate translate(tCoordinate in, uint dim) const {
		return ((in - bounds.coords[dim]) / bounds.dim[dim]) * (img.dim[dim] - img.coords[dim]);
	}

private:

	const dBox & bounds;
	const dBox & img;

	Cairo::RefPtr<Cairo::ImageSurface> cs;
	Cairo::RefPtr<Cairo::Context> cr;

};
