#pragma once
#include "IntersectionChecker.h"
#include <CGAL/Kernel_d/Point_d.h>
#include <CGAL/Kernel_d/Hyperplane_d.h>
#include "VectorOperations.h"

namespace LoadBalancing {

	template <uint D, typename Precision>
	struct HyperplaneIntersectionChecker : IntersectionChecker<D, Precision>
	{
		HyperplaneIntersectionChecker(dBox<D, Precision> bounds, dVector<D, Precision> normal,
		                              Precision offset)
			: IntersectionChecker<D, Precision>(std::move(bounds)), mNormal(std::move(normal)),
		      mOffset(offset) {}
		
		virtual bool intersects(const dSphere<D, Precision>& sphere) const
		{
			return distanceToPlane(sphere.center, mNormal, mOffset) > 0;
		}
		
		virtual std::unique_ptr<IntersectionChecker<D, Precision>> copy() const
		{
			return std::make_unique<HyperplaneIntersectionChecker>(*this);
		}
		
	private:
		dVector<D, Precision> mNormal;
		Precision mOffset;
	};
}
