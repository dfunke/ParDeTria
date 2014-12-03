#pragma once

#include <cairo/cairo.h>
#include <cairomm/cairomm.h>

#include <algorithm>

#include "Geometry.h"
#include "Logger.h"

class Painter {

public:
	typedef std::tuple<tCoordinate, tCoordinate, tCoordinate> tRGB;

public:

	Painter() {
		//non working version;
	}

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
	void draw(const dPoints & points);

	void draw(const dSimplex & simplex, const dPoints & points, bool drawInfinite = false);
	void draw(const dSimplices & simplices, const dPoints & points, bool drawInfinite = false);

	void drawCircumCircle(const dSimplex & simplex, const dPoints & points, bool drawInfinite = false);
	void drawCircumCircle(const dSimplices & simplices, const dPoints & points, bool drawInfinite = false);

	void drawPartition(const dPoints & points);

	void drawNeighbors(const dSimplex & simplex, const dSimplices & neighbors, const dPoints & points, bool drawInfinite = false);
	void drawNeighbors(const dSimplices & simplices, const dSimplices & neighbors, const dPoints & points, bool drawInfinite = false);

	void setColor(const tRGB & rgb, tCoordinate alpha = 1.0){
		setColor(std::get<0>(rgb), std::get<1>(rgb), std::get<2>(rgb), alpha);
	}

	void setColor(tCoordinate r, tCoordinate g, tCoordinate b, tCoordinate alpha = 1.0){
		cr->set_source_rgba(r, g, b, alpha);
	}

	void setLineWidth(tCoordinate w){
		cr->set_line_width(w);
	}

	void savePNG(const std::string & file) const {
		cs->write_to_png(file);
	}

public:

	static tRGB tetradicColor(uint i) {
		return tRGB(!(i & 2), i & 2, i & 1);
	}

private:

	void _init(const dBox & _bounds, uint _resolution);
	void _copy(const Painter & a);

	void line(const dPoints & points, uint a, uint b, bool drawInfinite);

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

#include <boost/progress.hpp>
#include <map>

class PainterBulkWriter {

public:

	~PainterBulkWriter(){

		if(!m_items.empty()){
			LOG << "Writing " << m_items.size() << " items" << std::endl;

			boost::progress_display progress(m_items.size(), LOG);

			for(auto & item : m_items){
				item.second.savePNG(item.first);
				++progress;
			}
		}
	}

	Painter & operator[](const std::string & file){
		return m_items[file];
	}

private:
	std::map<std::string, Painter> m_items;

};
