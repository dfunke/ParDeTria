#pragma once
#include "Geometry.h"

namespace LoadBalancing
{
	template <uint D, typename Precision>
	struct BoundsChecker : IntersectionChecker<D, Precision>
	{
		BoundsChecker(dBox<D, Precision> bounds)
			: IntersectionChecker(std::move(bounds)) {
		}

		virtual bool intersects(const dSphere<D, Precision>& sphere)
		{
			return bounds.intersects(sphere);
		}
	};
}
