#pragma once
#include "Geometry.h"
#include "load_balancing/IntersectionChecker.h"

namespace LoadBalancing
{
	template <uint D, typename Precision>
	struct BoundsIntersectionChecker : IntersectionChecker<D, Precision>
	{
		BoundsIntersectionChecker(dBox<D, Precision> bounds)
			: IntersectionChecker<D, Precision>(std::move(bounds)) {
		}

		virtual bool intersects(const dSphere<D, Precision>& sphere) const override
		{
			return IntersectionChecker<D, Precision>::bounds().intersects(sphere);
		}
	};
}
