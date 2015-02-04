#include "Partitioner.h"

template <uint D>
Partitioning<D>
dPartitioner<D>::partition(const Ids &ids, const dPoints<D> &points,
                           __attribute((unused))
                           const std::string &provenance) const {
  // do mid-point based partitioning for now
  auto stats = getPointStats(ids.begin(), ids.end(), points);

  PLOG << "Midpoint is " << stats.mid << std::endl;

  Partitioning<D> partitioning;
  for (uint i = 0; i < pow(2, D); ++i) {
    partitioning[i].id = i;

    for (uint d = 0; d < D; ++d) {
      partitioning[i].bounds.low[d] =
          i & (1 << d) ? stats.mid.coords[d] : stats.min.coords[d];
      partitioning[i].bounds.high[d] =
          i & (1 << d) ? stats.max.coords[d] : stats.mid.coords[d];
    }
  }

#ifndef NDEBUG
  auto inPartition = [&](const dPoint<D> &p, const uint partition) -> bool {
    for (uint i = 0; i < D; ++i)
      if (!(partitioning[partition].bounds.low[i] <= p.coords[i] &&
            p.coords[i] <= partitioning[partition].bounds.high[i]))
        return false;
    return true;
  };
#endif

  for (auto &id : ids) {
    ASSERT(points.contains(id));
    const auto &p = points[id];

    if (!p.isFinite())
      continue; // skip infinite points, they will be added later

    uint part = 0;
    for (uint dim = 0; dim < D; ++dim) {
      part |= (p.coords[dim] > stats.mid.coords[dim]) << dim;
    }

    ASSERT(inPartition(p, part));

    // LOG << "Adding " << p << " to " << part << std::endl;
    partitioning[part].points.insert(p.id);
  }

  // add infinite points
  for (uint i = 0; i < partitioning.size(); ++i) {
    for (uint k = dPoint<D>::cINF; k != 0; ++k) {
      partitioning[i].points.insert(k);
    }
  }

  return partitioning;
}

template <uint D>
Partitioning<D>
kPartitioner<D>::partition(const Ids &ids, const dPoints<D> &points,
                           __attribute((unused))
                           const std::string &provenance) const {
  // do mid-point based partitioning for now
  auto stats = getPointStats(ids.begin(), ids.end(), points);

  PLOG << "Midpoint is " << stats.mid << std::endl;

  Partitioning<D> partitioning;

  partitioning[0].id = 0;
  for (uint d = 0; d < D; ++d) {
    partitioning[0].bounds.low[d] = stats.min.coords[d];
    partitioning[0].bounds.high[d] =
        d == k ? stats.mid.coords[d] : stats.max.coords[d];
  }

  partitioning[1].id = 1;
  for (uint d = 0; d < D; ++d) {
    partitioning[1].bounds.low[d] =
        d == k ? stats.mid.coords[d] : stats.min.coords[d];
    partitioning[1].bounds.high[d] = stats.max.coords[d];
  }

#ifndef NDEBUG
  auto inPartition = [&](const dPoint<D> &p, const uint partition) -> bool {
    for (uint i = 0; i < D; ++i)
      if (!(partitioning[partition].bounds.low[i] <= p.coords[i] &&
            p.coords[i] <= partitioning[partition].bounds.high[i]))
        return false;
    return true;
  };
#endif

  for (auto &id : ids) {
    ASSERT(points.contains(id));
    const auto &p = points[id];

    if (!p.isFinite())
      continue; // skip infinite points, they will be added later

    uint part = (p.coords[k] > stats.mid.coords[k]);

    ASSERT(inPartition(p, part));

    // LOG << "Adding " << p << " to " << part << std::endl;
    partitioning[part].points.insert(p.id);
  }

  // add infinite points
  for (uint i = 0; i < partitioning.size(); ++i) {
    for (uint k = dPoint<D>::cINF; k != 0; ++k) {
      partitioning[i].points.insert(k);
    }
  }

  return partitioning;
}

template <uint D>
Partitioning<D>
CyclePartitioner<D>::partition(const Ids &ids, const dPoints<D> &points,
                               const std::string &provenance) const {
  // do mid-point based partitioning for now
  auto stats = getPointStats(ids.begin(), ids.end(), points);

  // cycle is lenght of provenance - 1 modulo D
  uint k = (provenance.size() - 1) % D;

  PLOG << "Midpoint is " << stats.mid << std::endl;
  PLOG << "Splitting dimension is " << k << std::endl;

  Partitioning<D> partitioning;

  partitioning[0].id = 0;
  for (uint d = 0; d < D; ++d) {
    partitioning[0].bounds.low[d] = stats.min.coords[d];
    partitioning[0].bounds.high[d] =
        d == k ? stats.mid.coords[d] : stats.max.coords[d];
  }

  partitioning[1].id = 1;
  for (uint d = 0; d < D; ++d) {
    partitioning[1].bounds.low[d] =
        d == k ? stats.mid.coords[d] : stats.min.coords[d];
    partitioning[1].bounds.high[d] = stats.max.coords[d];
  }

#ifndef NDEBUG
  auto inPartition = [&](const dPoint<D> &p, const uint partition) -> bool {
    for (uint i = 0; i < D; ++i)
      if (!(partitioning[partition].bounds.low[i] <= p.coords[i] &&
            p.coords[i] <= partitioning[partition].bounds.high[i])) {

        return false;
      }
    return true;
  };
#endif

  for (auto &id : ids) {
    ASSERT(points.contains(id));
    const auto &p = points[id];

    if (!p.isFinite())
      continue; // skip infinite points, they will be added later

    uint part = (p.coords[k] > stats.mid.coords[k]);

    ASSERT(inPartition(p, part));

    // LOG << "Adding " << p << " to " << part << std::endl;
    partitioning[part].points.insert(p.id);
  }

  // add infinite points
  for (uint i = 0; i < partitioning.size(); ++i) {
    for (uint k = dPoint<D>::cINF; k != 0; ++k) {
      partitioning[i].points.insert(k);
    }
  }

  return partitioning;
}

// specializations
template class dPartitioner<2>;
template class kPartitioner<2>;
template class CyclePartitioner<2>;

template class dPartitioner<3>;
template class kPartitioner<3>;
template class CyclePartitioner<3>;
