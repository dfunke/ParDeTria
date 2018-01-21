#pragma once
#include "Geometry.h"

namespace LoadBalancing
{
	template <uint D, typename Precision>
	struct IntersectionChecker
	{
		IntersectionChecker(dBox<D, Precision> bounds) : mBounds(bounds) {}
		virtual bool intersects(const dSphere<D, Precision>& sphere) const = 0;
		const dBox<D, Precision>& bounds() const { return mBounds; }
	private:
		dBox<D, Precision> mBounds;
	};
}	
