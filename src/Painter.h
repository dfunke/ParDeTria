#pragma once

#include <cairo/cairo.h>
#include <cairomm/cairomm.h>

#include <algorithm>

#include "Geometry.h"
#include "Logger.h"

class Painter {

public:

	Painter(const dBox & _bounds, uint _resolution = 10) {

		_init(_bounds, _resolution);
	}

	Painter(const Painter & a){

		_copy(a);
	}

	Painter & operator=(const Painter & a) {

		_copy(a);

		return *this;

	}

	void draw(const dPoint & point);
	void draw(const dPoints & points, bool drawPartition = true);
	void draw(const dSimplices & simplices, const dPoints & points, bool drawCC = false);

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

	void _init(const dBox & _bounds, uint _resolution);
	void _copy(const Painter & a);

	const dPoint & getPoint(const dPoints & points, uint a) const;

	void line(const dPoints & points, uint a, uint b);

	void circumCircle(const dSimplex & s, const dPoints & points);

	inline tCoordinate width() const {
		return imgDim(0);
	}

	inline tCoordinate height() const {
		return imgDim(1);
	}

	inline tCoordinate imgDim(const uint dim) const {
		return img.dim[dim] - img.coords[dim];
	}

	inline tCoordinate translatePoint(tCoordinate in, uint dim) const {
		return ((offset.coords[dim] + in - bounds.coords[dim]) / (3 * bounds.dim[dim])) * imgDim(dim);
	}

	inline tCoordinate translateLength(tCoordinate in, uint dim) const {
		return (in / (3 * bounds.dim[dim])) * imgDim(dim);
	}

private:

	dBox  bounds;
	dBox  img;
	dPoint offset;

	Cairo::RefPtr<Cairo::ImageSurface> cs;
	Cairo::RefPtr<Cairo::Context> cr;

};
