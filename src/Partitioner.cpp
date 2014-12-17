#include "Partitioner.h"

std::string to_string(const Partition &p) {
  std::stringstream o;
  o << p.id << " " << p.bounds << " [";
  for (auto it = p.points.begin(); it != p.points.end(); ++it) {
    if (it != p.points.begin())
      o << ", ";
    o << *it;
  }
  o << "]";
  return o.str();
}

std::ostream &operator<<(std::ostream &o, const Partition &p) {
  return o << to_string(p);
}

Partitioning dPartitioner::partition(const Ids &ids, const dPoints &points,
                                     __attribute((unused))
                                     const std::string &provenance) const {
  // do mid-point based partitioning for now
  auto stats = getPointStats(ids.begin(), ids.end(), points);

  PLOG << "Midpoint is " << stats.mid << std::endl;

  Partitioning partitioning;
  for (uint i = 0; i < pow(2, D); ++i) {
    partitioning[i].id = i;

    for (uint d = 0; d < D; ++d) {
      partitioning[i].bounds.coords[d] =
          i & (1 << d) ? stats.mid.coords[d] : stats.min.coords[d];
      partitioning[i].bounds.dim[d] =
          i & (1 << d) ? stats.max.coords[d] - stats.mid.coords[d]
                       : stats.mid.coords[d] - stats.min.coords[d];
    }
  }

#ifndef NDEBUG
  auto inPartition = [&](const dPoint &p, const uint partition) -> bool {
    for (uint i = 0; i < D; ++i)
      if (!(partitioning[partition].bounds.coords[i] <= p.coords[i] &&
            p.coords[i] <= partitioning[partition].bounds.coords[i] +
                               partitioning[partition].bounds.dim[i]))
        return false;
    return true;
  };
#endif

  for (auto &id : ids) {
    assert(points.contains(id));
    const auto &p = points[id];

    if (!p.isFinite())
      continue; // skip infinite points, they will be added later

    uint part = 0;
    for (uint dim = 0; dim < D; ++dim) {
      part |= (p.coords[dim] > stats.mid.coords[dim]) << dim;
    }

    assert(inPartition(p, part));

    // LOG << "Adding " << p << " to " << part << std::endl;
    partitioning[part].points.insert(p.id);
  }

  // add infinite points
  for (uint i = 0; i < partitioning.size(); ++i) {
    for (uint k = dPoint::cINF; k != 0; ++k) {
      partitioning[i].points.insert(k);
    }
  }

  return partitioning;
}

Partitioning kPartitioner::partition(const Ids &ids, const dPoints &points,
                                     __attribute((unused))
                                     const std::string &provenance) const {
  // do mid-point based partitioning for now
  auto stats = getPointStats(ids.begin(), ids.end(), points);

  PLOG << "Midpoint is " << stats.mid << std::endl;

  Partitioning partitioning;

  partitioning[0].id = 0;
  for (uint d = 0; d < D; ++d) {
    partitioning[0].bounds.coords[d] = stats.min.coords[d];
    partitioning[0].bounds.dim[d] =
        d == k ? stats.mid.coords[d] - stats.min.coords[d]
               : stats.max.coords[d] - stats.min.coords[d];
  }

  partitioning[1].id = 0;
  for (uint d = 0; d < D; ++d) {
    partitioning[1].bounds.coords[d] =
        d == k ? stats.mid.coords[d] : stats.min.coords[d];
    partitioning[1].bounds.dim[d] =
        d == k ? stats.max.coords[d] - stats.mid.coords[d]
               : stats.max.coords[d] - stats.min.coords[d];
  }

#ifndef NDEBUG
  auto inPartition = [&](const dPoint &p, const uint partition) -> bool {
    for (uint i = 0; i < D; ++i)
      if (!(partitioning[partition].bounds.coords[i] <= p.coords[i] &&
            p.coords[i] <= partitioning[partition].bounds.coords[i] +
                               partitioning[partition].bounds.dim[i]))
        return false;
    return true;
  };
#endif

  for (auto &id : ids) {
    assert(points.contains(id));
    const auto &p = points[id];

    if (!p.isFinite())
      continue; // skip infinite points, they will be added later

    uint part = (p.coords[k] > stats.mid.coords[k]);

    assert(inPartition(p, part));

    // LOG << "Adding " << p << " to " << part << std::endl;
    partitioning[part].points.insert(p.id);
  }

  // add infinite points
  for (uint i = 0; i < partitioning.size(); ++i) {
    for (uint k = dPoint::cINF; k != 0; ++k) {
      partitioning[i].points.insert(k);
    }
  }

  return partitioning;
}

Partitioning CyclePartitioner::partition(const Ids &ids, const dPoints &points,
                                         const std::string &provenance) const {
  // do mid-point based partitioning for now
  auto stats = getPointStats(ids.begin(), ids.end(), points);

  // cycle is lenght of provenance - 1 modulo D
  uint k = (provenance.size() - 1) % D;

  PLOG << "Midpoint is " << stats.mid << std::endl;
  PLOG << "Splitting dimension is " << k << std::endl;

  Partitioning partitioning;

  partitioning[0].id = 0;
  for (uint d = 0; d < D; ++d) {
    partitioning[0].bounds.coords[d] = stats.min.coords[d];
    partitioning[0].bounds.dim[d] =
        d == k ? stats.mid.coords[d] - stats.min.coords[d]
               : stats.max.coords[d] - stats.min.coords[d];
  }

  partitioning[1].id = 0;
  for (uint d = 0; d < D; ++d) {
    partitioning[1].bounds.coords[d] =
        d == k ? stats.mid.coords[d] : stats.min.coords[d];
    partitioning[1].bounds.dim[d] =
        d == k ? stats.max.coords[d] - stats.mid.coords[d]
               : stats.max.coords[d] - stats.min.coords[d];
  }

#ifndef NDEBUG
  auto inPartition = [&](const dPoint &p, const uint partition) -> bool {
    for (uint i = 0; i < D; ++i)
      if (!(partitioning[partition].bounds.coords[i] <= p.coords[i] &&
            p.coords[i] <= partitioning[partition].bounds.coords[i] +
                               partitioning[partition].bounds.dim[i]))
        return false;
    return true;
  };
#endif

  for (auto &id : ids) {
    assert(points.contains(id));
    const auto &p = points[id];

    if (!p.isFinite())
      continue; // skip infinite points, they will be added later

    uint part = (p.coords[k] > stats.mid.coords[k]);

    assert(inPartition(p, part));

    // LOG << "Adding " << p << " to " << part << std::endl;
    partitioning[part].points.insert(p.id);
  }

  // add infinite points
  for (uint i = 0; i < partitioning.size(); ++i) {
    for (uint k = dPoint::cINF; k != 0; ++k) {
      partitioning[i].points.insert(k);
    }
  }

  return partitioning;
}
