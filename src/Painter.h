#pragma once

#include <cairo/cairo.h>

#include "Geometry.h"

class Painter {

public:

	Painter(const dBox & _bounds, const dBox & _img)
	: bounds(_bounds), img(_img) {

		cs = cairo_image_surface_create(CAIRO_FORMAT_ARGB32,
				width(), height());
		cr = cairo_create(cs);

		//draw background white
		cairo_set_source_rgb(cr, 1,1,1);
		cairo_paint(cr);

		cairo_set_line_width(cr, 1.0);
		cairo_set_source_rgb(cr, 0,0,0);
	}

	Painter(const Painter & a)
	: bounds(a.bounds), img(a.img) {

		cs = cairo_image_surface_create(CAIRO_FORMAT_ARGB32,
				width(), height());
		cr = cairo_create(cs);

		//copy a's surface
		cairo_set_source_surface (cr, a.cs, 0, 0);
		cairo_paint (cr);
	}

	~Painter() {
		cairo_surface_destroy(cs);
	}

	void draw(const dPoints & points, bool drawPartition = true){

		for(const auto & p : points){
			cairo_arc(cr, translate(p.coords[0], 0), translate(p.coords[1],1), 5, 0, 2*M_PI);
			cairo_fill(cr);
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
			cairo_move_to(cr, translate(midpoint.coords[0], 0),0);
			cairo_line_to(cr, translate(midpoint.coords[0], 0),height());

			cairo_move_to(cr, 0, translate(midpoint.coords[1],1));
			cairo_line_to(cr, width(), translate(midpoint.coords[1],1));

			cairo_stroke(cr);
		}
	}

	void draw(const dSimplices & simplices, const dPoints & points){

		for(const auto & s : simplices){
			for(uint d = 0; d < D; ++d){
				line(points, s.vertices[d], s.vertices[d+1]);
			}
			//close the loop
			line(points, s.vertices[D], s.vertices[0]);
		}

		cairo_stroke(cr);
	}

	void setColor(float r, float g, float b){
		if(r > 1 || g > 1 || b > 1){
			//rescale colors
			r /= 255; g /= 255; b /= 255;
		}

		cairo_set_source_rgb(cr, r, g, b);
	}

	void setLineWidth(float w){
		cairo_set_line_width(cr, w);
	}

	void savePNG(const std::string & file) const {
		cairo_surface_write_to_png(cs, file.c_str());
	}

private:

	void line(const dPoints & points, uint a, uint b){
		if(a != INF && b  != INF){
			cairo_move_to(cr, translate(points[a].coords[0],0), translate(points[a].coords[1],1));
			cairo_line_to(cr, translate(points[b].coords[0],0), translate(points[b].coords[1],1));
		}
	}

	float width() const {
		return img.dim[0] - img.coords[0];
	}

	float height() const {
		return img.dim[1] - img.coords[1];
	}

	float translate(float in, uint dim) const {
		return ((in - bounds.coords[dim]) / bounds.dim[dim]) * (img.dim[dim] - img.coords[dim]);
	}

private:

	const dBox & bounds;
	const dBox & img;

	cairo_surface_t * cs;
	cairo_t * cr;

};
