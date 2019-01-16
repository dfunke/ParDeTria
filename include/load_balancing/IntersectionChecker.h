#pragma once
#include <functional>
#include "Geometry.h"

namespace LoadBalancing
{
	template <uint D, typename Precision>
	struct IntersectionChecker
	{
		IntersectionChecker(dBox<D, Precision> bounds) : mBounds(bounds) {}
		virtual ~t IntersectionChecker() = default;
		virtual bool intersects(const dSphere<D, Precision>& sphere) const = 0;
		virtual std::unique_ptr<IntersectionChecker<D, Precision>> copy() const = 0;
		const dBox<D, Precision>& bounds() const { return mBounds; }
	private:
		dBox<D, Precision> mBounds;
	};
	
	template <uint D, typename Precision>
	struct IntersectionPartition {
		Point_Ids pointIds;
		std::unique_ptr<IntersectionChecker<D, Precision>> intersectionChecker;
	};

	using PartitionAssigner = std::function<size_t(tIdType)>;

	template <uint D, typename Precision>
	using IntersectionPartitionMakerFunction =
		std::function<std::vector<IntersectionPartition<D, Precision>>(const dPoints<D, Precision>&,
																	   const Point_Ids&,
																	   size_t,
																	   PartitionAssigner)>;

	template <uint D, typename Precision>
	struct IntersectionPartitionMaker
	{
	};
}	
